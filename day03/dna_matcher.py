import argparse
from dataclasses import dataclass
from pathlib import Path
import re

MINIMUM_MATCH_PERCENT = 80.0
DNA_PATTERN = re.compile(r"^[ACGT]+$")


@dataclass(frozen=True)
class Difference:
    position: int
    expected: str
    actual: str


@dataclass(frozen=True)
class Match:
    start: int
    end: int
    candidate: str
    match_percent: float
    differences: list[Difference]


def validate_dna(sequence: str, description: str) -> str:
    # הסרת כל הרווחים, הטאבים וירידות השורה בתוך הרצף
    clean_sequence = re.sub(r"\s+", "", sequence).upper()
    if not clean_sequence or not DNA_PATTERN.match(clean_sequence):
        raise ValueError(f"{description} must contain only A, C, G, and T")
    return clean_sequence


def find_matches(sequence: str, query: str) -> list[Match]:
    clean_seq = validate_dna(sequence, "The DNA sequence")
    clean_query = validate_dna(query, "The query sequence")

    q_len = len(clean_query)
    s_len = len(clean_seq)
    if q_len > s_len:
        return []

    max_allowed_mismatches = int(q_len * (1 - MINIMUM_MATCH_PERCENT / 100))
    matches: list[Match] = []

    # Sliding window פשוט ומהיר
    for i in range(s_len - q_len + 1):
        window = clean_seq[i : i + q_len]

        differences: list[Difference] = []
        mismatch_count = 0

        for offset, (expected, actual) in enumerate(zip(clean_query, window)):
            if expected != actual:
                mismatch_count += 1
                if mismatch_count > max_allowed_mismatches:
                    break
                differences.append(Difference(i + offset + 1, expected, actual))
        else:
            match_percent = (q_len - mismatch_count) / q_len * 100
            matches.append(
                Match(
                    start=i + 1,
                    end=i + q_len,
                    candidate=window,
                    match_percent=match_percent,
                    differences=differences,
                )
            )

    return sorted(matches, key=lambda m: (-m.match_percent, m.start))


def print_matches(source_name: str, sequence: str, query: str) -> None:
    matches = find_matches(sequence, query)
    print(f"\n{source_name}:")
    if not matches:
        print(f"No candidates matched at least {MINIMUM_MATCH_PERCENT:g}%")
        return

    for index, match in enumerate(matches):
        label = "best match" if index == 0 else "match"
        print(
            f"{label}: {match.match_percent:g}% match. starts at position: {match.start} "
            f"and ends at position: {match.end}. candidate: {match.candidate}"
        )
        if match.differences:
            changes = ", ".join(
                f"position {diff.position}: {diff.expected} -> {diff.actual}"
                for diff in match.differences
            )
            print(f"  changes: {changes}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find DNA subsequences similar to a query sequence."
    )
    parser.add_argument(
        "query",
        nargs="?",
        help="DNA sequence to search for (candidates are the same length)",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        type=Path,
        help="DNA text file to search (prompted for if omitted)",
    )
    args = parser.parse_args()

    try:
        input_file = args.input_file or Path(input("DNA file path: ").strip())
        query = args.query or input("Wanted DNA sequence: ")

        sequence = input_file.read_text(encoding="utf-8")
        print_matches(input_file.name, sequence, query)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    user_input = input("DO you want to try another Sequence or DNA? (press Enter to exit press Y to continue): ")
    user_input = user_input.strip().lower()
    if user_input == "Y" or user_input == "y":
        main()


if __name__ == "__main__":
    main()