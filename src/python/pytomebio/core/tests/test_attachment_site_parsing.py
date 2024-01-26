import pytest

from pytomebio.core.attachment_site import AttachmentSite


@pytest.mark.parametrize(
    "string_to_parse,expected_attachment_site",
    [
        (
            "test_site_palindrome:AAA:GT:TCGA",
            AttachmentSite("test_site_palindrome", "AAA", "GT", "TCGA"),
        ),
        (
            "test_site_complement:AAA:GT:AAA",
            AttachmentSite("test_site_complement", "AAA", "GT", "TTT"),
        ),
        (
            "test_site_reverse_complement:AAA:GT:CTA",
            AttachmentSite("test_site_reverse_complement", "AAA", "GT", "TAG"),
        ),
    ],
)
def test_right_branch(string_to_parse: str, expected_attachment_site: AttachmentSite) -> None:
    assert AttachmentSite.parse_attachment_site(site=string_to_parse) == expected_attachment_site


@pytest.mark.parametrize(
    "input_strings_that_raise_assertion_errors",
    [
        ("name:AAA:GT:AAA:AAA"),
        ("name:AAA:GT:AAA-AAA"),
        ("name:Alice::Bob"),
        (""),
        ("name"),
        ("name:AAA"),
        ("name:AAA:GT"),
        ("name-AAA-GT-AAA"),
        ("name,AAA,GT,AAA"),
    ],
)
def test_greater_than_three_colons_in_input_string(
    input_strings_that_raise_assertion_errors: str,
) -> None:
    with pytest.raises(AssertionError):
        AttachmentSite.parse_attachment_site(site=input_strings_that_raise_assertion_errors)
