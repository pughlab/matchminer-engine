"""Copyright 2016 Dana-Farber Cancer Institute"""

__author__ = 'priti,james,zachary'
import os
import networkx as nx

from matchengine.engine import MatchEngine
from matchengine.utilities import build_oncotree
from tests import TestSetUp

YAML_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/yaml/'))


def read_file(file):
    fh = open(file, 'r')
    data = ''
    for stream in fh:
        data += stream
    return data


class TestMatchEngine(TestSetUp):

    # global vars.
    me = None

    def setUp(self):
        super(TestMatchEngine, self).setUp()

        # add clinical collection
        self.db.clinical.insert_many([{
            'ONCOTREE_PRIMARY_DIAGNOSIS_NAME': onc,
            'SAMPLE_ID': self.sample_id,
            'VITAL_STATUS': 'alive',
            'MRN': self.mrn,
            'REPORT_DATE': self.static_date,
            'BIRTH_DATE': self.static_date,
            'GENDER': 'Male'
        } for onc in ['_LIQUID_', '_SOLID_']])

        self.add_clinical()
        self.add_genomic()

    def tearDown(self):
        self.db.clinical.drop()
        self.db.genomic.drop()

    def test_validate_yaml_format(self):

        # Read yaml file
        test_inp = read_file(os.path.join(YAML_DIR, '00-000.yml'))

        # Parse yaml document, The status should be 1 as validation should fail
        status, data = self.me.validate_yaml_format(test_inp)
        assert status == 1

        # Parse yaml document, The status should be 0 as validation should pass
        test_inp = read_file(os.path.join(YAML_DIR, '00-001.yml'))
        status, data = self.me.validate_yaml_format(test_inp)
        assert status == 0

        # Assert the protocol id is correctly annotated
        assert data['protocol_no'] == '00-001'
        assert data['protocol_id'] == 00001

    def test_validate_yaml_data(self):

        # Read yaml file
        test_inp = read_file(os.path.join(YAML_DIR, '00-001.yml'))
        status, data = self.me.validate_yaml_format(test_inp)

        # Validate the schema for yaml document
        errors = self.me.validate_yaml_data(data)

        # Assert that there are no errors
        assert len(errors) == 0

        test_inp = read_file(os.path.join(YAML_DIR, '00-002.yml'))
        status, data = self.me.validate_yaml_format(test_inp)

        # Assert Schema check fails
        errors = self.me.validate_yaml_data(data)
        assert errors['protocol_id'][0] == 'required field'

        # assert we don't need a match clause at root.
        test_inp = read_file(os.path.join(YAML_DIR, '00-003.yml'))
        status, data = self.me.validate_yaml_format(test_inp)
        assert status == 0

    def test_run_query(self):

        # reinstantiate MatchEngine so that the set of all sample ids in the database includes the documents that were
        # posted by the unit test setUp
        self.me = MatchEngine(self.db)

        # define a clinical node
        node = {'type': 'clinical', 'value': {'ONCOTREE_PRIMARY_DIAGNOSIS': 'Adrenal Gland', 'AGE_NUMERICAL': '>=18'}}

        # run query
        result, matches = self.me.run_query(node)

        # assert number of samples are 89
        assert len(result) == 1, len(result)
        assert self.sample_id in result

        # define a genomic criteria
        node = {'type': 'genomic', 'value': {'HUGO_SYMBOL': '!BRAF'}}

        # run query
        result, matches = self.me.run_query(node)

        # assert number of samples is 1
        assert len(result) == 9, len(result)
        assert self.sample_id not in result

        self.add_genomic_v2()
        node = {'type': 'genomic', 'value': {'HUGO_SYMBOL': 'WHSC1'}}
        result, matches = self.me.run_query(node)
        assert 'actionability' in matches[0]
        assert matches[0]['mmr_status'] == 'Proficient (MMR-P / MSS)'

    def test_prepare_clinical_criteria(self):

        onc = 'ONCOTREE_PRIMARY_DIAGNOSIS'
        oncname = 'ONCOTREE_PRIMARY_DIAGNOSIS_NAME'

        # create clinical criteria
        item = {'AGE_NUMERICAL': '>=18', onc: 'Melanoma'}

        # convert to mongo query
        c = self.me.prepare_clinical_criteria(item)

        # assert length of ONCOTREE_DIAGNOSIS for 'Melanoma' is 9
        assert len(c[oncname]['$in']) == 8

        # check 'BIRTH DATE mongo query does less than search'
        assert c['BIRTH_DATE'].keys()[0] == '$lte'

        # check "!"
        item = {'AGE_NUMERICAL': '<=18', onc: '!Melanoma'}
        c = self.me.prepare_clinical_criteria(item)
        assert c['BIRTH_DATE'].keys()[0] == '$gte'
        assert len(c[oncname]['$nin']) == 8

        # check with a list of diagnoses
        item = {'AGE_NUMERICAL': '>=18', onc: ['!Melanoma', '!Glioblastoma', 'Pheochromocytoma', 'Astrocytoma']}
        c = self.me.prepare_clinical_criteria(item)
        assert c['BIRTH_DATE'].keys()[0] == '$lte'
        assert '$nin' in c[oncname]
        assert '$in' in c[oncname]
        assert len(c[oncname]['$in']) == 2
        assert len(c[oncname]['$nin']) == 12

        # check _SOLID_ && _LIQUID_
        liquid_item = {'AGE_NUMERICAL': '>=18', onc: '_LIQUID_'}
        solid_item = {'AGE_NUMERICAL': '>=18', onc: '_SOLID_'}
        liqc = self.me.prepare_clinical_criteria(liquid_item)
        solc = self.me.prepare_clinical_criteria(solid_item)
        assert len(liqc[oncname]['$in']) == 51, len(liqc[oncname]['$in'])
        assert len(solc[oncname]['$in']) == 561, len(solc[oncname]['$in'])

    def test_prepare_genomic_criteria(self):

        # create a genomic criteria
        item = {'HUGO_SYMBOL': '!KRAS', 'PROTEIN_CHANGE': 'p.V600E'}

        # convert to mongo query
        c, neg, _ = self.me.prepare_genomic_criteria(item)

        # check ! symbol
        assert c['$and'][0]['TRUE_HUGO_SYMBOL']['$eq'] == 'KRAS'
        assert neg is True

        # check protein change
        assert c['$and'][0]['TRUE_PROTEIN_CHANGE']['$eq'] == 'p.V600E'

        # check wildtype
        assert '$or' in c['$and'][1], c
        assert c['$and'][1]['$or'] == [{'WILDTYPE': False}, {'WILDTYPE': {'$exists': False}}]

    def test_sv(self):

        # create a genomic criteria
        item = {'HUGO_SYMBOL': 'KRAS', 'VARIANT_CATEGORY': 'SV'}

        # convert to mongo query
        c, neg, _ = self.me.prepare_genomic_criteria(item)

        # check sv
        assert c['$and'][0]['VARIANT_CATEGORY']['$eq'] == 'SV'

        # check wildtype
        assert '$or' in c['$and'][1], c
        assert c['$and'][1]['$or'] == [{'WILDTYPE': False}, {'WILDTYPE': {'$exists': False}}]

    def test_create_trial_tree(self):
        # parse yaml file and create tree.
        test_inp = read_file(os.path.join(YAML_DIR, '00-004.yml'))
        status, trial_tree = self.me.create_trial_tree(test_inp)

        # assert it was successful.
        assert status == 0

        # check that we have 2 match trees.
        cnt = 0
        for n in trial_tree.nodes():
            if 'match_tree' in trial_tree.node[n]:
                cnt += 1
        assert cnt == 2

        # parse yaml file and create tree.
        test_inp = read_file(os.path.join(YAML_DIR, '00-004.yml'))
        status, trial_tree = self.me.create_trial_tree(test_inp)

        # check that we have 14 nodes.
        assert len(list(trial_tree.nodes())) == 4

        # assert we have 3 match trees.
        cnt = 0
        for n in trial_tree.nodes():
            if 'match_tree' in trial_tree.node[n]:
                cnt += 1
        assert cnt == 2

    def test_create_match_tree(self):

        # parse yaml file and create trial tree.
        test_inp = read_file(os.path.join(YAML_DIR, '00-004.yml'))
        status, trial_tree = self.me.create_trial_tree(test_inp)

        # get the graph for the only match.
        g = None
        i = -1
        number_of_nodes = [2, 10]
        edges = [[(1, 2)], [(1, 2), (2, 3), (2, 4), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (4, 10)]]
        for n in trial_tree.nodes():
            if 'match_tree' in trial_tree.node[n]:
                i += 1
                g = trial_tree.node[n]['match_tree']
                assert g is not None

                # Check if tree contain correct number of nodes
                assert g.number_of_nodes() == number_of_nodes[i]

                # Check edges are correctly created
                assert list(nx.dfs_edges(g)) == edges[i]

                # First node should be and
                if i == 0:
                    assert g.node[1]['type'] == 'match'
                    assert g.node[2]['type'] == 'clinical'
                    print g.node[2]['value']
                    assert g.node[2]['value'] == {'disease_status': ['Advanced']}

                elif i == 1:
                    assert g.node[1]['type'] == 'match'
                    assert g.node[2]['type'] == 'and'
                    assert g.node[3]['type'] == 'genomic'
                    assert g.node[3]['value'] == {'hugo_symbol': 'IDH1', 'wildcard_protein_change': 'p.R132', 'variant_category': 'Mutation'}
                    assert g.node[4]['type'] == 'or'
                    assert g.node[5]['type'] == 'clinical'
                    assert g.node[5]['value'] == {'age_numerical': '>=18', 'oncotree_primary_diagnosis': '_SOLID_'}
                    assert g.node[6]['type'] == 'clinical'
                    assert g.node[6]['value'] == {'age_numerical': '>=18', 'oncotree_primary_diagnosis': 'Diffuse Glioma'}
                    assert g.node[7]['type'] == 'clinical'
                    assert g.node[7]['value'] == {'age_numerical': '>=18', 'oncotree_primary_diagnosis': 'Encapsulated Glioma'}
                    assert g.node[8]['type'] == 'clinical'
                    assert g.node[8]['value'] == {'age_numerical': '>=18', 'oncotree_primary_diagnosis': 'Cholangiocarcinoma'}
                    assert g.node[9]['type'] == 'clinical'
                    assert g.node[9]['value'] == {'age_numerical': '>=18', 'oncotree_primary_diagnosis': 'Acute Myeloid Leukemia'}
                    assert g.node[10]['type'] == 'clinical'
                    assert g.node[10]['value'] == {'age_numerical': '>=18', 'oncotree_primary_diagnosis': 'Myelodysplasia'}

    def test_search_oncotree_diagnosis(self):

        # Glioblastoma
        c = {'ONCOTREE_PRIMARY_DIAGNOSIS_NAME': {'$eq': 'Glioblastoma'}}
        oncotree = build_oncotree()
        c['ONCOTREE_PRIMARY_DIAGNOSIS_NAME'] = self.me._search_oncotree_diagnosis(oncotree, c)
        conc = c['ONCOTREE_PRIMARY_DIAGNOSIS_NAME']
        assert '$in' in conc, conc
        assert conc['$in'] == ['Small Cell Glioblastoma', 'Gliosarcoma', 'Glioblastoma Multiforme', 'Glioblastoma'], conc

        # Melanoma
        c = {'ONCOTREE_PRIMARY_DIAGNOSIS_NAME': {'$eq': 'Melanoma'}}
        conc = self.me._search_oncotree_diagnosis(oncotree, c)
        assert '$in' in conc, conc
        assert conc['$in'] == [
            'Melanoma', 'Congenital Nevus', 'Genitourinary Mucosal Melanoma', 'Cutaneous Melanoma',
            'Melanoma of Unknown Primary', 'Desmoplastic Melanoma', 'Lentigo Maligna Melanoma', 'Acral Melanoma'
        ]
