"""A module for translating between alignments and edits sequences."""


def get_edits(p: str, q: str) -> tuple[str, str, str]:
    """Extract the edit operations from a pairwise alignment.

    Args:
        p (str): The first row in the pairwise alignment.
        q (str): The second row in the pairwise alignment.

    Returns:
        str: The list of edit operations as a string.

    >>> get_edits('ACCACAGT-CATA', 'A-CAGAGTACAAA')
    ('ACCACAGTCATA', 'ACAGAGTACAAA', 'MDMMMMMMIMMMM')

    """
    assert len(p) == len(q)
    # FIXME: do the actual calculations here
    orig, trans, cigar = "","",""

    for i in range(len(p)):
        if p[i] != "-" and q[i] != "-":
            orig += f"{p[i]}"
            trans += f"{q[i]}"
            cigar += f"M"
        elif p[i] == "-" and q[i] != "-":
            trans += f"{q[i]}"
            cigar += f"I"
        elif p[i] != "-" and q[i] == "-":
            orig += f"{p[i]}"
            cigar += f"D"

    return orig, trans, cigar


def local_align(p: str, x: str, i: int, edits: str) -> tuple[str, str]:
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> local_align("ACCACAGTCATA", "GTACAGAGTACAAA", 2, "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')

    """
    # FIXME: Compute the alignment rows
    return align(p, x[i:], edits)


def align(p: str, q: str, edits: str) -> tuple[str, str]:
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The first sequence to align.
        q (str): The second sequence to align
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> align("ACCACAGTCATA", "ACAGAGTACAAA", "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')

    """
    # FIXME: Compute the alignment rows
    orig = ""
    trans = ""
    pCounter, qCounter = 0 , 0
    for i in edits:
        if i == "M":
            orig += f"{p[pCounter]}"; pCounter += 1
            trans += f"{q[qCounter]}"; qCounter += 1
        elif i == "D":
            orig += f"{p[pCounter]}"; pCounter += 1
            trans += f"-"
        elif i == "I":
            orig += f"-"
            trans += f"{q[qCounter]}"; qCounter += 1


    return orig, trans


def edit_dist(p: str, x: str, i: int, edits: str) -> int:
    """Get the distance between p and the string that starts at x[i:]
    using the edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        int: The distance from p to x[i:?] described by edits

    >>> edit_dist("accaaagta", "cgacaaatgtcca", 2, "MDMMIMMMMIIM")
    5
    """
    # FIXME: Compute the edit distance
    distance = 0
    firstRow, secondRow = local_align(p, x, i, edits)
    firstCounter, secondCounter = 0,0
    for k in range(len(edits)):
        if edits[k] == "D" or edits[k] == "I":
            distance += 1
        else:
            if not firstRow[k] == secondRow[k]:
                distance += 1

    return distance
