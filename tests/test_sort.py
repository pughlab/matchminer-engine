"""Copyright 2016 Dana-Farber Cancer Institute"""

from matchengine.sort import *
from tests import TestSetUp


class TestSort(TestSetUp):

    def setUp(self):
        super(TestSort, self).setUp()
        pass

    def tearDown(self):
        pass

    def test_add_sort_value(self):

        li = []
        li = add_sort_value(8, 0, li)
        assert li == [8]

        li = add_sort_value(0, 0, li)
        assert li == [0]

        li = add_sort_value(1, 1, li)
        assert li == [0, 1]

        li = [0, 1, 0, 0]
        li = add_sort_value(1, 0, li)
        assert li == [0, 1, 0, 0]

    def test_sort_by_tier(self):

        sort_order = {('01', 'p01'): []}
        match = {
            'mrn': '01',
            'sample_id': '01',
            'protocol_no': 'p01',
            'mmr_status': 'MMR-Deficient',
            'tier': None,
            'actionability': None,
            'wildtype': False,
            'variant_category': 'MUTATION'
        }
        sort_order = sort_by_tier(match, sort_order)
        assert sort_order[('01', 'p01')][0] == 0

        sort_order = {('01', 'p01'): []}
        match['mmr_status'] = None
        match['tier'] = 1
        sort_order = sort_by_tier(match, sort_order)
        assert sort_order[('01', 'p01')][0] == 1

        sort_order = {('01', 'p01'): []}
        match['tier'] = 2
        sort_order = sort_by_tier(match, sort_order)
        assert sort_order[('01', 'p01')][0] == 2

        sort_order = {('01', 'p01'): []}
        match['tier'] = 3
        sort_order = sort_by_tier(match, sort_order)
        assert sort_order[('01', 'p01')][0] == 4

        sort_order = {('01', 'p01'): []}
        match['tier'] = 4
        sort_order = sort_by_tier(match, sort_order)
        assert sort_order[('01', 'p01')][0] == 5

        sort_order = {('01', 'p01'): []}
        match['tier'] = None
        sort_order = sort_by_tier(match, sort_order)
        assert sort_order[('01', 'p01')][0] == 7

        sort_order = {('01', 'p01'): []}
        match['variant_category'] = 'CNV'
        sort_order = sort_by_tier(match, sort_order)
        assert sort_order[('01', 'p01')][0] == 3

        sort_order = {('01', 'p01'): []}
        match['tier'] = None
        match['variant_category'] = None
        match['wildtype'] = None
        match['mmr_status'] = None
        match['clinical_only'] = True
        sort_order = sort_by_tier(match, sort_order)
        assert sort_order[('01', 'p01')][0] == 7

    def test_sort_by_match_type(self):

        sort_order = {('01', 'p01'): [0]}
        match = {
            'mrn': '01',
            'sample_id': '01',
            'protocol_no': 'p01',
            'match_type': None
        }
        sort_order = sort_by_match_type(match, sort_order)
        assert sort_order[('01', 'p01')][1] == 2

        match['match_type'] = 'gene'
        sort_order = sort_by_match_type(match, sort_order)
        assert sort_order[('01', 'p01')][1] == 1

        match['match_type'] = 'variant'
        sort_order = sort_by_match_type(match, sort_order)
        assert sort_order[('01', 'p01')][1] == 0

    def test_sort_by_cancer_type(self):

        sort_order = {('01', 'p01'): [0, 0]}
        match = {
            'mrn': '01',
            'sample_id': '01',
            'protocol_no': 'p01',
            'cancer_type_match': 'all_solid'
        }
        sort_order = sort_by_cancer_type(match, sort_order)
        assert sort_order[('01', 'p01')][2] == 1

        match['cancer_type_match'] = 'all_liquid'
        sort_order = sort_by_cancer_type(match, sort_order)
        assert sort_order[('01', 'p01')][2] == 1

        match['cancer_type_match'] = 'specific'
        sort_order = sort_by_cancer_type(match, sort_order)
        assert sort_order[('01', 'p01')][2] == 0

        match['cancer_type_match'] = 'unknown'
        sort_order = sort_by_cancer_type(match, sort_order)
        assert sort_order[('01', 'p01')][2] == 0

    def test_sort_by_coordinating_center(self):

        sort_order = {('01', 'p01'): [0, 0, 0]}
        match = {
            'mrn': '01',
            'sample_id': '01',
            'protocol_no': 'p01',
            'coordinating_center': 'Massachussetts General Hospital'
        }
        sort_order = sort_by_coordinating_center(match, sort_order)
        assert sort_order[('01', 'p01')][3] == 1

        match['coordinating_center'] = 'Dana-Farber Cancer Institute'
        sort_order = sort_by_coordinating_center(match, sort_order)
        assert sort_order[('01', 'p01')][3] == 0

    def test_sort_by_reverse_protocol_no(self):

        sort_order = {
            ('01', '11-111'): [0, 0, 0, 0],
            ('01', '09-999'): [0, 0, 0, 1],
            ('01', '15-000'): [7, 0, 0, 0],
            ('01', '15-111'): [7, 0, 0, 0],
            ('01', '22-222'): [7, 1, 0, 0],
        }
        matches = [
            {'protocol_no': '11-111', 'sample_id': '01'},
            {'protocol_no': '09-999', 'sample_id': '01'},
            {'protocol_no': '15-000', 'sample_id': '01'},
            {'protocol_no': '15-111', 'sample_id': '01'},
            {'protocol_no': '22-222', 'sample_id': '01'},
            {'protocol_no': '22-222', 'sample_id': '01'}
        ]
        sort_order = sort_by_reverse_protocol_no(matches, sort_order)
        assert sort_order[('01', '11-111')] == [0, 0, 0, 0, 3]
        assert sort_order[('01', '09-999')] == [0, 0, 0, 1, 4]
        assert sort_order[('01', '15-000')] == [7, 0, 0, 0, 2]
        assert sort_order[('01', '15-111')] == [7, 0, 0, 0, 1]
        assert sort_order[('01', '15-111')] == [7, 0, 0, 0, 1]
        assert sort_order[('01', '15-111')] == [7, 0, 0, 0, 1]
        assert sort_order[('01', '22-222')] == [7, 1, 0, 0, 0]

    def test_final_sort(self):

        mso = {}
        sort_order = {
            ('01', '11-111'): [0, 0, 0, 0, 1],
            ('01', '09-999'): [0, 0, 0, 1, 0],
            ('01', '12-000'): [2, 3, 1, 0, 2],
            ('01', '12-222'): [2, 3, 1, 1, 3],
            ('01', '12-333'): [2, 3, 1, 0, 4],
            ('01', '13-000'): [2, 0, 0, 0, 5],
            ('01', '15-333'): [2, 3, 1, 0, 6]
        }
        mso = final_sort(sort_order, mso)
        assert mso[('01', '11-111')] == 0
        assert mso[('01', '09-999')] == 1
        assert mso[('01', '13-000')] == 2
        assert mso[('01', '12-000')] == 3
        assert mso[('01', '12-333')] == 4
        assert mso[('01', '15-333')] == 5
        assert mso[('01', '12-222')] == 6

    def test_add_sort_order(self):

        tm = self.get_demo_trial_matches()
        tm = add_sort_order(tm)

        # print tm[['protocol_no', 'sort_order']].sort_values(by='sort_order', ascending=True)
        assert tm[['protocol_no', 'sort_order']].sort_values(by='sort_order', ascending=True).protocol_no.tolist() == \
            [
                '0003-000',  # tm13 (SV match (gets a sort order of -1))
                '0001-000',  # tm11 (mmr status deficient)
                '111-000',  # tm1  (tier 1, variant match, specific cancer type, DFCI, higher protocol #)
                '000-000',   # tm10 (tier 1, variant match, specific cancer type, DFCI, lower protocol #)
                '999-000',   # tm9  (tier 1, variant match, specific cancer type, MGH)
                '888-000',   # tm8  (tier 1, variant match, solid cancer type)
                '777-000',   # tm7  (tier 1, gene match)
                '444-000',   # tm4  (tier 2)
                '333-000',   # tm3  (CNV)
                '555-000',   # tm5  (tier 3)
                '666-000',   # tm6  (tier 4, higher protocol #)
                '222-000',   # tm2  (tier 4)
                '0002-000',  # tm12 (wildtype)
                '0004-000',  # tm14 (clinical only)
            ]
