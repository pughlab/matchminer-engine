"""Copyright 2016 Dana-Farber Cancer Institute"""

import os
import json
import yaml
import random
import string
import unittest
import datetime as dt
from bson.objectid import ObjectId

from matchengine.utilities import get_db
from matchengine.engine import MatchEngine

YAML_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/yaml/'))


class TestSetUp(unittest.TestCase):

    def setUp(self):
        """
        Descriptions of test patients

        1: >18, Adrenal Gland, Female, BRAF F346R Mutation
        2: >18, Melanoma, Female, EGFR L858R Mutation
        3: >18, Melanoma, Female, EGFR F346A Mutation
        4: >18, Melanoma, Female, EGFR F346B Mutation
        5: >18, Melanoma, Female, EGFR F000F Mutation
        6: >0.5 && <18, Melanoma, Male, EGFR SV
        7: >0.5 && <18, Glioblastoma, Male, EGFR CNV Hetero del
        8: >0.5 && <18, Glioblastoma, Male, EGFR CNV Gain
        9: >0.5 && <18, Glioblastoma, Male, EGFR CNV Homo del
        10: <0.5, Glioblastoma, Male, EGFR CNV High amp

        Descriptions of test trials
        00-001.yml: dose: EGFR L858R && >=18/_SOLID_
        00-002.yml: arm: EGFR L858R && >=18/_SOLID_
        00-003.yml: step: EGFR L858R && >=18/_SOLID_
        00-004.yml dose: EGFR L858R && >=18/_SOLID_
        00-005.yml 2 doses: EGFR L858R && >=18/_SOLID_
        00-006.yml exon: !13
        """

        self.db = get_db(None)
        for res in ["clinical", "dashboard", "filter", "genomic", "hipaa", "match", "normalize", "oplog"
                    "response", "statistics", "status", "team", "trial", "trial_match", "user"]:
            self.db.drop_collection(res)

        self.me = MatchEngine(self.db)

        self.trials = {}
        self.clinical_id = ObjectId()
        self.mrn = 'TCGA-BH-A1FR'
        self.sample_id = 'TCGA-OR-A5J1'
        self.mrns = [self.mrn] + [self.__random_id() for _ in range(9)]
        self.sample_ids = [self.sample_id] + [self.__random_id() for _ in range(9)]
        self.clinical_ids = [self.clinical_id] + [ObjectId() for _ in range(9)]
        self.static_date = dt.datetime.today()

        # clinical collection
        self.oncotree_diagnoses = ['Adrenal Gland'] + ['Melanoma'] * 5 + ['Glioblastoma'] * 4
        self.genders = ['Female'] * 5 + ['Male'] * 5

        # ages
        adult = self.static_date - dt.timedelta(days=365*19)
        child = self.static_date - dt.timedelta(days=365*10)
        infant = self.static_date - dt.timedelta(days=30*4)
        self.ages = [adult] * 5 + [child] * 4 + [infant]

        self.clinical = [{
            '_id': clinical_id,
            'ONCOTREE_PRIMARY_DIAGNOSIS_NAME': diagnosis,
            'SAMPLE_ID': sample_id,
            'VITAL_STATUS': 'alive',
            'MRN': mrn,
            'REPORT_DATE': self.static_date,
            'BIRTH_DATE': age,
            'GENDER': gender
        } for diagnosis, gender, age, clinical_id, sample_id, mrn in zip(
            self.oncotree_diagnoses, self.genders, self.ages, self.clinical_ids, self.sample_ids, self.mrns)]

        # genomic collection
        self.genes = ['BRAF'] + ['EGFR'] * 9
        self.protein_changes = ['p.F346R', 'p.L858R', 'p.F346A', 'p.F346B', 'p.F000F', None, None, None, None, None]
        self.variant_categories = ['MUTATION'] * 5 + ['SV', 'CNV', 'CNV', 'CNV', 'CNV']
        self.wildtypes = [False] * 10
        self.cnv_calls = [None, None, None, None, None, None,
                          'Heterozygous deletion', 'Gain', 'Homozygous deletion', 'High level amplification']
        self.genomic = [{
            'TRUE_VARIANT_CLASSIFICATION': 'In_Frame_Del',
            'TRUE_PROTEIN_CHANGE': protein_change,
            'VARIANT_CATEGORY': variant_category,
            'CHROMOSOME': 'chr3',
            'POSITION': 178952085,
            'TRUE_STRAND': '+',
            'WILDTYPE': wildtype,
            'CLINICAL_ID': _id,
            'CNV_CALL': cnv_call,
            'TRUE_HUGO_SYMBOL': gene,
            'SAMPLE_ID': sample_id,
            'TRUE_TRANSCRIPT_EXON': 19
        } for protein_change, variant_category, wildtype, cnv_call, gene, _id, sample_id in zip(
            self.protein_changes, self.variant_categories, self.wildtypes,
            self.cnv_calls, self.genes, self.clinical_ids, self.sample_ids
        )]

        # test trials
        self.test_trials = ['00-001', '00-002', '00-003']

        # demo match results
        pnos = ['00-001', '00-001', '00-001', '00-002', '00-002', '00-002']
        mlevels = ['arm', 'arm', 'arm', 'dose', 'dose', 'dose']
        iids = ['1', '2', '3', '4', '5', '6']
        galts = ['Alt1', 'Alt2', 'Alt2', 'Alt3', 'Alt3', 'Alt3']
        self.matches = [{
            'mrn': 'SAMPLE1',
            'sample_id': 'SAMPLE1-ID',
            'protocol_no': protocol_no,
            'match_level': match_level,
            'internal_id': internal_id,
            'genomic_alteration': genomic_alteration
        } for protocol_no, match_level, internal_id, genomic_alteration in zip(
            pnos, mlevels, iids, galts
        )]

    def tearDown(self):
        self.db.map.drop()
        self.db.trial_match.drop()

    def add_clinical(self):
        """Add all clinical documents to database needed for unit tests"""
        self.db.clinical.insert_many(self.clinical)

    def add_genomic(self):
        """Add all genomic documents to database needed for unit tests"""
        self.db.genomic.insert_many(self.genomic)

    def add_genomic_v2(self):
        """Adds a genomic with OncoPanel layout version 2"""
        genomic = {
            'TRUE_VARIANT_CLASSIFICATION': None,
            'TRUE_PROTEIN_CHANGE': None,
            'VARIANT_CATEGORY': 'CNV',
            'CHROMOSOME': None,
            'POSITION': None,
            'TRUE_STRAND': None,
            'WILDTYPE': False,
            'CLINICAL_ID': ObjectId(),
            'CNV_CALL': 'Heterozygous deletion',
            'TRUE_HUGO_SYMBOL': 'WHSC1',
            'SAMPLE_ID': 'FAKE-01',
            'TRUE_TRANSCRIPT_EXON': None,
            'ACTIONABILITY': 'actionable',
            'MMR_STATUS': 'Proficient (MMR-P / MSS)'
        }
        self.db.genomic.insert(genomic)

    def add_genomic_for_exon_mutation(self):
        """
        One of the unit tests has to test for the scenario where a trial calls for no mutation in a specific exon
        This method sets up the db for that unit test
        """

        # reset
        self.db.genomic.drop()

        # add genomic collection to db
        variants = ['MUTATION', 'MUTATION', 'CNV']
        clinical_ids = [ObjectId() for _ in range(3)]
        cnv_calls = [None, None, 'Homozygous deletion']
        sample_ids = ['1', '2', '3']                        # 2 and 3 should match; 1 should not match
        exons = [13, 19, 13]
        g = [{
            'TRUE_VARIANT_CLASSIFICATION': 'In_Frame_Del',
            'TRUE_PROTEIN_CHANGE': 'p.V600E',
            'VARIANT_CATEGORY': variant_category,
            'CHROMOSOME': 'chr3',
            'POSITION': 178952085,
            'TRUE_STRAND': '+',
            'WILDTYPE': False,
            'CLINICAL_ID': _id,
            'CNV_CALL': cnv_call,
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'SAMPLE_ID': sample_id,
            'TRUE_TRANSCRIPT_EXON': exon
        } for variant_category, _id, cnv_call, sample_id, exon in zip(
            variants, clinical_ids, cnv_calls, sample_ids, exons
        )]

        # add
        self.db.genomic.insert_many(g)

    def add_genomic_for_regex(self):
        """Adds genomic entries to specifically test pymongo $regex operations"""

        # reset
        self.db.genomic.drop()

        # add genomic collection to db
        # p.A0
        muts = ['p.A000Z', 'p.B0_A0B', 'p.B0A', 'p.A0B', 'p.A0fs*6', 'p.A0*', 'p.A0_B12insL']
        g = [{'TRUE_PROTEIN_CHANGE': mut} for mut in muts]
        self.db.genomic.insert_many(g)
        return muts

    def add_trials(self, trials=None):
        """Add all trial documents to database needed for unit tests"""

        if not trials:
            trials=self.test_trials

        # read trials from direcotry
        for yml in os.listdir(YAML_DIR):
            ymlpath = os.path.join(YAML_DIR, yml)

            if yml.split('.')[0] not in trials:
                continue

            # convert yml to json format
            with open(ymlpath) as f:
                t = yaml.load(f.read())
                self.trials[t['protocol_no']] = t

                # add trial to db
                self.db.trial.insert_one(t)

    def add_wildtype(self):
        """Adds a variety of wildtype genomic entries to the database"""

        clinical_ids = [ObjectId() for _ in range(3)]
        wts = [True, True, False]
        sample_ids = ['1', '2', '3']
        g = [{
            'TRUE_VARIANT_CLASSIFICATION': 'In_Frame_Del',
            'TRUE_PROTEIN_CHANGE': 'p.V600E',
            'VARIANT_CATEGORY': 'MUTATION',
            'CHROMOSOME': 'chr3',
            'POSITION': 178952085,
            'TRUE_STRAND': '+',
            'WILDTYPE': wt,
            'CLINICAL_ID': _id,
            'CNV_CALL': None,
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'SAMPLE_ID': sample_id,
            'TRUE_TRANSCRIPT_EXON': 13
        } for _id, sample_id, wt in zip(clinical_ids, sample_ids, wts)]

        g.append({
            'TRUE_VARIANT_CLASSIFICATION': 'In_Frame_Del',
            'TRUE_PROTEIN_CHANGE': 'p.V600E',
            'VARIANT_CATEGORY': 'MUTATION',
            'CHROMOSOME': 'chr3',
            'POSITION': 178952085,
            'TRUE_STRAND': '+',
            'CLINICAL_ID': ObjectId(),
            'CNV_CALL': None,
            'TRUE_HUGO_SYMBOL': 'EGFR',
            'SAMPLE_ID': '4',
            'TRUE_TRANSCRIPT_EXON': 13
        })

        self.db.genomic.insert_many(g)

    def add_msi(self):

        clinical_ids = [ObjectId() for _ in range(3)]
        mmr_statuses = ['Proficient (MMR-P / MSS)', 'Deficient (MMR-D / MSI-H)', 'Indeterminate (see note)']
        sample_ids = ['1', '2', '3']
        g = [{
            'TRUE_VARIANT_CLASSIFICATION': None,
            'TRUE_PROTEIN_CHANGE': None,
            'VARIANT_CATEGORY': 'SIGNATURE',
            'CHROMOSOME': None,
            'POSITION': None,
            'TRUE_STRAND': None,
            'WILDTYPE': False,
            'CLINICAL_ID': _id,
            'CNV_CALL': None,
            'TRUE_HUGO_SYMBOL': None,
            'SAMPLE_ID': sample_id,
            'TRUE_TRANSCRIPT_EXON': 13,
            'MMR_STATUS': mmr
        } for _id, sample_id, mmr in zip(clinical_ids, sample_ids, mmr_statuses)]
        self.db.genomic.insert_many(g)

    @staticmethod
    def get_demo_trial_matches():

        trial_match = {
            "code": "1",
            "match_type": "variant",
            "vital_status": "alive",
            "true_cdna_change": "c.1799T>A",
            "sample_id": "111",
            "chromosome": "1",
            "oncotree_primary_diagnosis_name": "Cutaneous Melanoma",
            "internal_id": "01",
            "true_variant_classification": "Missense_Mutation",
            "ord_physician_name": "FAKE PHYSICIAN",
            "ord_physician_email": "fake_physician@fake.fake",
            "variant_category": "MUTATION",
            "protocol_no": "111-000",
            "wildtype": False,
            "canonical_strand": "-",
            "first_last": "FIRST LAST",
            "true_hugo_symbol": "BRAF",
            "trial_accrual_status": "open",
            "match_level": "step",
            "tier": 1,
            "allele_fraction": 0.38,
            "genomic_alteration": "BRAF p.V600E",
            "gender": "Female",
            "reference_allele": "A",
            "mrn": "111",
            "true_protein_change": "p.V600E",
            "true_transcript_exon": 15,
            "position": 140453136,
            "cnv_call": None,
            "actionability": None,
            "mmr_status": None,
            "cancer_type_match": "specific",
            "coordinating_center": "Dana-Farber Cancer Institute"
        }
        tm1 = trial_match.copy()
        tm2 = trial_match.copy()
        tm3 = trial_match.copy()
        tm4 = trial_match.copy()
        tm5 = trial_match.copy()
        tm6 = trial_match.copy()
        tm7 = trial_match.copy()
        tm8 = trial_match.copy()
        tm9 = trial_match.copy()
        tm10 = trial_match.copy()
        tm11 = trial_match.copy()
        tm12 = trial_match.copy()
        tm13 = trial_match.copy()
        tm14 = trial_match.copy()

        tm2['actionability'] = "actionable"
        tm2['tier'] = 4
        tm3['variant_category'] = "CNV"
        tm3['tier'] = None
        tm4['tier'] = 2
        tm5['tier'] = 3
        tm6['tier'] = 4
        tm7['match_type'] = 'gene'
        tm8['cancer_type_match'] = 'all_solid'
        tm9['coordinating_center'] = 'Massachussetts General Hospital'
        tm10['protocol_no'] = '11-111'
        tm11['mmr_status'] = 'Deficient (MMR-D / MSI-H)'
        tm12['wildtype'] = True
        tm12['tier'] = None
        tm13['genomic_alteration'] = ' Structural Variation'
        tm13['true_hugo_symbol'] = None
        tm14['tier'] = None
        tm14['mmr_status'] = None
        tm14['variant_category'] = None
        tm14['wildtype'] = None
        tm14['clinical_only'] = True

        tm2['protocol_no'] = '222-000'
        tm3['protocol_no'] = '333-000'
        tm4['protocol_no'] = '444-000'
        tm5['protocol_no'] = '555-000'
        tm6['protocol_no'] = '666-000'
        tm7['protocol_no'] = '777-000'
        tm8['protocol_no'] = '888-000'
        tm9['protocol_no'] = '999-000'
        tm10['protocol_no'] = '000-000'
        tm11['protocol_no'] = '0001-000'
        tm12['protocol_no'] = '0002-000'
        tm13['protocol_no'] = '0003-000'
        tm14['protocol_no'] = '0004-000'

        trial_matches = [tm1, tm2, tm3, tm4, tm5, tm6, tm7, tm8, tm9, tm10, tm11, tm12, tm13, tm14]
        return trial_matches

    @staticmethod
    def __random_id():
        """Returns a random MRN or SAMPLE_ID formatted "XXXX-XX-XXXX"""
        return 'TCGA-%s-%s' % (
            ''.join(random.choice(string.ascii_uppercase) for _ in range(2)),
            ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
        )

    @staticmethod
    def _debug(t):
        try:
            return json.dumps(t, sort_keys=True, indent=4)
        except TypeError:
            return t
