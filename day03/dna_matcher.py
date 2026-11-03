import argparse
import re
from pathlib import Path


MINIMUM_MATCH_PERCENT = 80
DNA_PATTERN = re.compile(r"[ACGT]+", re.IGNORECASE)


def validate_dna(sequence: str, description: str) -> str:
    sequence = sequence.strip()
    if not DNA_PATTERN.fullmatch(sequence):
        raise ValueError(f"{description} must contain only A, C, G, and T")
    return sequence.upper()


def find_matches(
    sequence: str, query: str
) -> list[tuple[int, str, float, list[tuple[int, str, str]]]]:
    sequence = validate_dna(sequence, "The DNA sequence")
    query = validate_dna(query, "The query sequence")
    if len(query) > len(sequence):
        return []

    candidates = re.finditer(rf"(?=([ACGT]{{{len(query)}}}))", sequence)
    matches = []
    for candidate in candidates:
        window = candidate.group(1)
        differences = [
            (candidate.start() + index + 1, expected, actual)
            for index, (expected, actual) in enumerate(zip(query, window))
            if expected != actual
        ]
        match_percent = (len(query) - len(differences)) / len(query) * 100
        if match_percent >= MINIMUM_MATCH_PERCENT:
            matches.append((candidate.start() + 1, window, match_percent, differences))

    return sorted(matches, key=lambda match: (-match[2], match[0]))


def print_matches(source_name: str, sequence: str, query: str) -> None:
    matches = find_matches(sequence, query)
    print(f"\n{source_name}:")
    if not matches:
        print(f"No candidates matched at least {MINIMUM_MATCH_PERCENT}%")
        return

    for index, (start, candidate, match_percent, differences) in enumerate(matches):
        label = "best match" if index == 0 else "match"
        end = start + len(query) - 1
        print(
            f"{label}: {match_percent:g}% match. starts at position: {start} "
            f"and ends at position: {end}. candidate: {candidate}"
        )
        if differences:
            changes = ", ".join(
                f"position {position}: {expected} -> {actual}"
                for position, expected, actual in differences
            )
            print(f"  changes: {changes}")


def run_tests(query: str) -> None:
    test_folder = Path(__file__).parent / "test"
    test_files = sorted(test_folder.glob("*.txt"))
    if not test_files:
        raise FileNotFoundError(f"No .txt files found in {test_folder}")

    for test_file in test_files:
        sequence = test_file.read_text(encoding="utf-8")
        print_matches(test_file.name, sequence, query)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find DNA subsequences similar to a query sequence."
    )
    parser.add_argument(
        "query",
        help="DNA sequence to search for (candidates are the same length)",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        type=Path,
        help="DNA text file to search (searches the test folder if omitted)",
    )
    args = parser.parse_args()

    try:
        if args.input_file is None:
            run_tests(args.query)
        else:
            sequence = args.input_file.read_text(encoding="utf-8")
            print_matches(args.input_file.name, sequence, args.query)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
