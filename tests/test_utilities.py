"""Copyright 2016 Dana-Farber Cancer Institute"""

import copy

from matchengine.utilities import *
from matchengine.settings import months
from matchengine.engine import MatchEngine as me
from tests import TestSetUp


YAML_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/yaml/'))


class TestUtilities(TestSetUp):

    def setUp(self):
        super(TestUtilities, self).setUp()

        self.today = dt.datetime.today()

        # get map from database
        self.mapping = list(self.db.map.find())

        # add clinical entries
        self.add_clinical()

    def tearDown(self):
        self.db.clinical.drop()
        self.db.trial.drop()

    def test_get_sampleids_from_mrns(self):
        mrn_map = samples_from_mrns(self.db, [self.mrn])
        assert mrn_map == {self.sample_id: self.mrn}, 'MRN MAP: %s\nMRN: %s\nSAMPLE_ID: %s' % (
            mrn_map, self.mrn, self.sample_id)

    def test_get_months(self):

        self.today = self.today.replace(year=2016, month=11, day=3)

        # The return of "get_months" provides the month to replace and the years to subtract from today's date
        months_, year = get_months(float(.5), self.today)
        assert months_ == 5, months_
        assert year == 0, year

        months_, year = get_months(float(.25), self.today)
        assert months_ == 8, months_
        assert year == 0, year

        months_, year = get_months(float(1.5), self.today)
        assert months_ == 5, months_
        assert year == 1, year

        months_, year = get_months(float(10.25), self.today)
        assert months_ == 8, months_
        assert year == 10, year

        # A special use case is if subtracting months from today's date causes the date to proceed into the previous
        # year. At this point, we set today to February and assert that the correct math is applied.
        today = self.today.replace(month=2)
        months_, year = get_months(float(.5), today)
        assert months_ == 8, months_
        assert year == -1, year

        months_, year = get_months(float(10.25), today)
        assert months_ == 11, months_
        assert year == -11, year

    def test_search_birth_date(self):

        c = {'BIRTH_DATE': {'$eq': '<=18'}}
        self._assert_age(search_birth_date(c)['$gte'], 18)

        c = {'BIRTH_DATE': {'$eq': '>=18'}}
        self._assert_age(search_birth_date(c)['$lte'], 18)

        c = {'BIRTH_DATE': {'$eq': '<18'}}
        self._assert_age(search_birth_date(c)['$gt'], 18)

        c = {'BIRTH_DATE': {'$eq': '>18'}}
        self._assert_age(search_birth_date(c)['$lt'], 18)

        c = {'BIRTH_DATE': {'$eq': '>.5'}}
        self._assert_age(search_birth_date(c)['$lt'], 0, 6)

        # c = {'BIRTH_DATE': {'$eq': '<=10.25'}}
        # self._assert_age(search_birth_date(c)['$gte'], 10, 3)

    def test_normalize_fields(self):
        assert normalize_fields(self.mapping, 'age_numerical')[0] == 'BIRTH_DATE'
        assert normalize_fields(self.mapping, 'exon')[0] == 'TRUE_TRANSCRIPT_EXON'
        assert normalize_fields(self.mapping, 'hugo_symbol')[0] == 'TRUE_HUGO_SYMBOL'
        assert normalize_fields(self.mapping, 'protein_change')[0] == 'TRUE_PROTEIN_CHANGE'
        assert normalize_fields(self.mapping, 'wildcard_protein_change')[0] == 'TRUE_PROTEIN_CHANGE'
        assert normalize_fields(self.mapping, 'variant_classification')[0] == 'TRUE_VARIANT_CLASSIFICATION'
        assert normalize_fields(self.mapping, 'variant_category')[0] == 'VARIANT_CATEGORY'
        assert normalize_fields(self.mapping, 'cnv_call')[0] == 'CNV_CALL'
        assert normalize_fields(self.mapping, 'wildtype')[0] == 'WILDTYPE'

    def test_normalize_values(self):
        assert normalize_values(self.mapping, 'variant_category', 'Mutation') == ('VARIANT_CATEGORY', 'MUTATION')
        assert normalize_values(self.mapping, 'variant_category', 'Copy Number Variation') == ('VARIANT_CATEGORY', 'CNV')
        assert normalize_values(self.mapping, 'variant_category', 'Structural Variation') == ('VARIANT_CATEGORY', 'SV')
        assert normalize_values(self.mapping, 'variant_category', '!Mutation') == ('VARIANT_CATEGORY', '!MUTATION')
        assert normalize_values(self.mapping, 'cnv_call', 'High Amplification') == ('CNV_CALL', 'High level amplification')
        assert normalize_values(self.mapping, 'cnv_call', 'Homozygous Deletion') == ('CNV_CALL', 'Homozygous deletion')
        assert normalize_values(self.mapping, 'cnv_call', 'Heterozygous Deletion') == ('CNV_CALL', 'Heterozygous deletion')
        assert normalize_values(self.mapping, 'wildtype', 'true') == ('WILDTYPE', True)
        assert normalize_values(self.mapping, 'wildtype', 'false') == ('WILDTYPE', False)

    def test_build_oncotree(self):
        onco_tree = build_oncotree()
        assert onco_tree.nodes()

    def test_build_gquery(self):
        # wildcard protein change
        key, txt, neg, _ = build_gquery('wildcard_protein_change', 'p.F346')
        assert txt == '^p.F346[A-Z]'
        assert key == '$regex'
        assert neg is False

        key, txt, neg, _ = build_gquery('wildcard_protein_change', 'F346')
        assert txt == '^p.F346[A-Z]'
        assert key == '$regex'
        assert neg is False

        key, txt, neg, _ = build_gquery('wildcard_protein_change', '!p.F346')
        assert txt == '^p.F346[A-Z]'
        assert key == '$regex'
        assert neg is True

        # not equal to
        key, txt, neg, _ = build_gquery('protein_change', '!p.F346')
        assert txt == 'p.F346'
        assert key == '$eq'
        assert neg is True

        # equal to
        key, txt, neg, _ = build_gquery('protein_change', 'p.F346')
        assert txt == 'p.F346'
        assert key == '$eq'
        assert neg is False

        # not exon
        key, txt, neg, _ = build_gquery('exon', '!13')
        assert txt == 13
        assert key == '$eq'
        assert neg is True

        # any variation
        key, txt, neg, _ = build_gquery('variant_category', 'Any Variation')
        assert txt == ['MUTATION', 'CNV']
        assert key == '$in'
        assert neg is False

        # mmr status
        key, txt, neg, sv = build_gquery('mmr_status', 'MMR-Proficient')
        assert txt == 'Proficient (MMR-P / MSS)'
        assert key == '$eq'
        assert neg is False
        assert sv is False

        key, txt, neg, sv = build_gquery('mmr_status', 'MMR-Deficient')
        assert txt == 'Deficient (MMR-D / MSI-H)'
        assert key == '$eq'
        assert neg is False
        assert sv is False

        # ms status
        key, txt, neg, sv = build_gquery('ms_status', 'MSI-H')
        assert txt == 'Deficient (MMR-D / MSI-H)'
        assert key == '$eq'
        assert neg is False
        assert sv is False

        key, txt, neg, sv = build_gquery('ms_status', 'MSI-L')
        assert txt == 'Proficient (MMR-P / MSS)'
        assert key == '$eq'
        assert neg is False
        assert sv is False

        key, txt, neg, sv = build_gquery('ms_status', 'MSS')
        assert txt == 'Proficient (MMR-P / MSS)'
        assert key == '$eq'
        assert neg is False
        assert sv is False


    def test_build_cquery(self):

        field = 'ONCOTREE_PRIMARY_DIAGNOSIS_NAME'
        txt = 'Melanoma'

        # OncoTree primary diagnosis
        c = build_cquery({}, field, txt)
        assert c == {field: {'$eq': txt}}, c

        # Not equal to
        nottxt = '!Melanoma'
        c = build_cquery({}, field, nottxt)
        assert c == {field: {'$ne': txt}}, c

        # lists
        txt = ['Melanoma', 'Glioblastoma']
        c = build_cquery({}, field, txt)
        assert c == {field: {'$in': txt}}, c

        # not in lists
        nottxt = ['!Melanoma', '!Glioblastoma']
        c = build_cquery({}, field, nottxt)
        assert c == {field: {'$nin': txt}}, c

        # both
        bothtxt = ['!Melanoma', '!Glioblastoma', 'Pheochromocytoma', 'Astrocytoma']
        c = build_cquery({}, field, bothtxt)
        assert c == {field: {'$nin': txt, '$in': bothtxt[2:]}}, c

        # non expanding oncotree primary diagnosis
        txt = 'Peritoneum'
        c = build_cquery({}, field, txt)
        assert c == {field: {'$eq': txt}}, c

    def test_format_genomic_alteration(self):

        # Gene only
        item = {'TRUE_HUGO_SYMBOL': 'EGFR'}
        gquery = {'TRUE_HUGO_SYMBOL': {'$eq': 'EGFR'}}
        g, is_variant = format_genomic_alteration(item, gquery)
        assert g == 'EGFR', g
        assert is_variant == 'gene'

        # Protein Change
        item = {
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'TRUE_PROTEIN_CHANGE': 'p.V600E'
        }
        gquery = {'TRUE_HUGO_SYMBOL': {'$eq': 'EGFR'}, 'TRUE_PROTEIN_CHANGE': {'$eq': 'p.V600E'}}
        g, is_variant = format_genomic_alteration(item, gquery)
        assert g == 'EGFR p.V600E', g
        assert is_variant == 'variant'

        # CNV
        item = {
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'CNV_CALL': 'High level amplification'
        }
        gquery = {'TRUE_HUGO_SYMBOL': {'$eq': 'EGFR'}, 'TRUE_PROTEIN_CHANGE': {'$ne': 'p.V600E'}}
        g, is_variant = format_genomic_alteration(item, gquery)
        assert g == 'EGFR High level amplification', g
        assert is_variant == 'variant'

        # SV
        item = {
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'VARIANT_CATEGORY': 'SV'
        }
        gquery = {'$and': [
            {'TRUE_HUGO_SYMBOL': {'$eq': 'EGFR'}, 'TRUE_PROTEIN_CHANGE': {'$ne': 'p.V600E'}},
            {'$or': [{'WILDTYPE': False}, {'WILDTYPE': {'$exists': False}}]}
        ]}
        g, is_variant = format_genomic_alteration(item, gquery)
        assert g == 'EGFR Structural Variation', g
        assert is_variant == 'variant'

    def test_format_not_match(self):

        # ! HUGO only, single
        gquery = {'TRUE_HUGO_SYMBOL': {'$eq': 'BRAF'}}
        alt, is_variant = format_not_match(gquery)
        assert alt == '!BRAF', alt
        assert is_variant == 'gene'

        # ! HUGO only, multiple
        gquery = {'TRUE_HUGO_SYMBOL': {'$in': ['BRAF', 'EGFR']}}
        alt, is_variant = format_not_match(gquery)
        assert alt == '!BRAF, !EGFR', alt
        assert is_variant == 'gene'

        # ! HUGO with PROTEIN CHANGE
        gquery = {'TRUE_HUGO_SYMBOL': {'$eq': 'EGFR'}, 'TRUE_PROTEIN_CHANGE': {'$eq': 'p.V600E'}}
        alt, is_variant = format_not_match(gquery)
        assert alt == '!EGFR p.V600E', alt
        assert is_variant == 'variant'

        # ! HUGO with multiple PROTEIN CHANGES
        gquery = {'TRUE_HUGO_SYMBOL': {'$eq': 'EGFR'}, 'TRUE_PROTEIN_CHANGE': {'$in': ['p.V600E', 'p.V600B']}}
        alt, is_variant = format_not_match(gquery)
        assert alt == '!EGFR p.V600E, p.V600B', alt
        assert is_variant == 'variant'

        # ! multiple HUGOS with single PROTEIN CHANGE
        gquery = {'TRUE_HUGO_SYMBOL': {'$in': ['BRAF', 'EGFR']}, 'TRUE_PROTEIN_CHANGE': {'$eq': 'p.V600E'}}
        alt, is_variant = format_not_match(gquery)
        assert alt == '!BRAF, !EGFR p.V600E', alt
        assert is_variant == 'variant'

        # ! multiple HUGOS with multiple PROTEIN CHANGES
        gquery = {
            'TRUE_HUGO_SYMBOL': {'$in': ['BRAF', 'EGFR']},
            'TRUE_PROTEIN_CHANGE': {'$in': ['p.V600E', 'p.V600B']}
        }
        alt, is_variant = format_not_match(gquery)
        assert alt == '!BRAF, !EGFR p.V600E, p.V600B', alt
        assert is_variant == 'variant'

        # ! PROTEIN CHANGE regex
        gquery = {'TRUE_PROTEIN_CHANGE': {'$regex': '^p.V600[A-Z]'}}
        alt, is_variant = format_not_match(gquery)
        assert alt == '!p.V600', alt
        assert is_variant == 'variant'

        # ! HUGO with regex PROTEIN CHANGE
        gquery = {'TRUE_PROTEIN_CHANGE': {'$regex': '^p.V600[A-Z]'}, 'TRUE_HUGO_SYMBOL': {'$eq': 'BRAF'}}
        alt, is_variant = format_not_match(gquery)
        assert alt == '!BRAF p.V600', alt
        assert is_variant == 'variant'

    def test_format_query(self):

        query = {'$eq': 'p.V600E'}
        alteration = format_query(query)
        assert alteration == 'p.V600E', alteration

        query = {'$regex': '^p.V600[A-Z]'}
        alteration = format_query(query)
        assert alteration == 'p.V600', alteration

        query = {'$in': ['MUTATION', 'CNV', 'SV']}
        alteration = format_query(query)
        assert alteration == 'MUTATION, CNV, SV', alteration

        query = {'$eq': 'BRAF'}
        alteration = format_query(query, gene=True)
        assert alteration == '!BRAF', alteration

        query = {'$in': ['BRAF', 'KRAS']}
        alteration = format_query(query, gene=True)
        assert alteration == '!BRAF, !KRAS', alteration

    def test_get_cancer_type_match(self):

        trial = {}
        ctm = get_cancer_type_match(trial)
        assert ctm == 'unknown'

        trial = {'_summary': {'tumor_types': ['_SOLID_']}}
        ctm = get_cancer_type_match(trial)
        assert ctm == 'all_solid'

        trial = {'_summary': {'tumor_types': ['_LIQUID_']}}
        ctm = get_cancer_type_match(trial)
        assert ctm == 'all_liquid'

        trial = {'_summary': {'tumor_types': ['Lung']}}
        ctm = get_cancer_type_match(trial)
        assert ctm == 'specific'

    def test_get_coordinating_center(self):

        trial = {}
        tcc = get_coordinating_center(trial)
        assert tcc == 'unknown'

        trial = {'_summary': {'coordinating_center': 'Massachusetts General Hospital'}}
        tcc = get_coordinating_center(trial)
        assert tcc == 'Massachusetts General Hospital'

    def _assert_age(self, bd, age, month=None):

        if month:
            try:
                assert bd.strftime('%B') == months[self.today.month - month - 1]
            except AssertionError:
                assert bd.strftime('%B') == months[self.today.month - month]
            assert self.today.strftime('%d') == bd.strftime('%d')
        else:
            assert self.today.strftime('%B %d') == bd.strftime('%B %d')

        # Unit test will fail the assertion if it crosses the year boundary
        try:
            assert int(self.today.strftime('%Y')) == int(bd.strftime('%Y')) + age
        except AssertionError:
            assert int(self.today.strftime('%Y')) == int(bd.strftime('%Y')) + age + 1
