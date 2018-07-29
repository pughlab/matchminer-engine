import datetime as dt

from tests import TestSetUp


class TestSV(TestSetUp):

    def setUp(self):
        super(TestSV, self).setUp()

        # add genomic collection
        self.db.genomic.drop()
        self.db.clinical.drop()
        self.db.trial.drop()

        self.db.clinical.insert_one({
            "SAMPLE_ID": "MATCH",
            "MRN": "MRN00",
            "VITAL_STATUS": "alive",
            "BIRTH_DATE": dt.datetime.today().replace(year=1900, month=11, day=3),
            "ONCOTREE_PRIMARY_DIAGNOSIS_NAME": "Glioblastoma Multiforme"
        })

        self.db.genomic.insert_one({
            "SAMPLE_ID": "MATCH",
            "VARIANT_CATEGORY": "SV",
            "STRUCTURAL_VARIANT_COMMENT": "An ETV6-NTRK3 fusion is identified (chr12:12035285 to chr15:88559895). "
        })

        self.db.genomic.insert_one({
            "SAMPLE_ID": "MATCH",
            "VARIANT_CATEGORY": "SV",
            "STRUCTURAL_VARIANT_COMMENT": "An ETV6-BRAF fusion is identified (chr12:12035285 to chr15:88559895). "
        })

        self.db.trial.insert_one({
            "protocol_no": "00-000",
            "treatment_list": {
                "step": [
                    {
                        "step_internal_id": 000,
                        "step_code": "0",
                        "step_type": "0",
                        "arm": [
                            {
                                "arm_code": "0",
                                "dose_level": [],
                                "arm_description": "0",
                                "arm_suspended": "N",
                                "arm_internal_id": 0,
                                "match": [
                                    {
                                        "and": [
                                            {
                                                "or": [
                                                    {
                                                        "genomic": {
                                                            "hugo_symbol": "MET",
                                                            "exon": 14,
                                                            "variant_category": "Mutation"
                                                        }
                                                    },
                                                    {
                                                        "genomic": {
                                                            "hugo_symbol": "MET",
                                                            "exon": 14,
                                                            "variant_category": "Structural Variation"
                                                        }
                                                    }
                                                ]
                                            },
                                            {
                                                "clinical": {
                                                    "age_numerical": ">=18",
                                                    "oncotree_primary_diagnosis": "Non-Small Cell Lung Cancer"
                                                }
                                            }
                                        ]
                                    }
                                ]
                            },
                            {
                                "arm_code": "0",
                                "dose_level": [],
                                "arm_description": "0",
                                "arm_suspended": "N",
                                "arm_internal_id": 0,
                                "match": [
                                    {
                                        "and": [
                                            {
                                                "or": [
                                                    {
                                                        "genomic": {
                                                            "hugo_symbol": "NTRK1",
                                                            "variant_category": "Structural Variation"
                                                        }
                                                    },
                                                    {
                                                        "genomic": {
                                                            "hugo_symbol": "NTRK2",
                                                            "variant_category": "Structural Variation"
                                                        }
                                                    },
                                                    {
                                                        "genomic": {
                                                            "hugo_symbol": "NTRK3",
                                                            "variant_category": "Structural Variation"
                                                        }
                                                    }
                                                ]
                                            },
                                            {
                                                "clinical": {
                                                    "age_numerical": ">=18",
                                                    "oncotree_primary_diagnosis": "_SOLID_"
                                                }
                                            }
                                        ]
                                    }
                                ]
                            }
                        ]
                    }
                ]
            }
        })

    def tearDown(self):
        self.db.clinical.drop()
        self.db.genomic.drop()
        self.db.trial.drop()

    def test_sv(self):

        mrn_map = dict(zip(["MATCH"], ["MRN00"]))
        trial_matches = []
        trial = self.db.trial.find_one({'protocol_no': '00-000'})
        trial_segment = trial['treatment_list']['step'][0]['arm'][1]
        match_segment = 'arm'

        # add sample id to trial_matches dictionary
        t = self.me._assess_match(mrn_map, trial_matches, trial, trial_segment, match_segment, 'open')
        assert len(t) == 1
