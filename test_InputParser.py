from chemkin import *
def test_not_implement():
    try:
        reactions = chemkin.from_xml('test_xml/rxns_bad_not_implement.xml')
    except NotImplementedError as e:
        assert type(e) == NotImplementedError
        print(e)
    try:
        reactions = chemkin.from_xml('test_xml/rxns_bad_not_implement_reverse.xml')
    except NotImplementedError as e:
        assert type(e) == NotImplementedError
        print(e)

def test_bad_reactant():
    try:
        reactions = chemkin.from_xml('test_xml/rxns_bad_reactants.xml')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_correct():
    input_ = InputParser('test_xml/rxns.xml')
    assert len(input_) == 3
    assert len(eval(repr(input_))) == 3
    print(repr(input_))