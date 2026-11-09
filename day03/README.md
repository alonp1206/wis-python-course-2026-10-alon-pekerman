# DNA Sequence Matcher

`dna_matcher.py` searches a DNA sequence for sections that are similar to a
query sequence supplied by the user. Each candidate section has the same
length as the query. The program uses a regular expression lookahead to find
every possible starting position, including overlapping candidates.

## Requirements

- Python 3
- A query containing only DNA bases `A`, `C`, `G`, and `T`
- A text file containing one DNA sequence, or `.txt` files in the sibling
  `test` folder

The program accepts upper- or lowercase bases. It strips whitespace at the
beginning and end of a sequence and rejects any other characters.

## Run

From the repository's top-level folder, provide the query and the file path:

```powershell
python day03\dna_matcher.py ACGT day03\test\dna.txt
```

If either value is omitted, the program prompts for it. Running without
arguments prompts for the file path first, then the wanted DNA sequence:

```powershell
python day03\dna_matcher.py
```

The `run_tests(query)` function is also available for searching every `.txt`
file in the sibling `test` folder.

## Results

Only candidates with at least 80% exact base matches are shown. Results are
ordered by match percentage, best first, then by their starting position.
Positions are one-based and the end position is inclusive. For candidates
that are not exact matches, the output lists each differing position and the
change from the query base to the candidate base.

For example, searching for `ACGTA` could report:

```text
best match: 80% match. starts at position: 1 and ends at position: 5. candidate: ACGTT
  changes: position 5: A -> T
```

With a four-base query, a candidate must match all four bases to meet the 80%
threshold. Longer queries can show partial matches.

## Prompts used

### first

> write the following python3 program under folder "day03":
> mission: finds sub strings from a dna sequance which match a given sequance (given by the user when run)
> 1. the function needs to read a txt file
> - the file will contain only a string of (dna) letters without anything else
> 2. the length of the dna sequance can vary
> 3. it uses regular expressions for it
> 4. it should show a flow of canidates sequances which will show match of at least 80%
> - the output should look like:
> best match: x% of match. starts at position: **(position number)** and end at position: **(position number)**
> - and on
> - for less then 100% match is also needs to indicate where it varies (position and change)
> 5. it should also contain a function named: "run_tests" which will read all the txt files from the "test" folder (will always sits in the same folder of the program).

### second

> add a README.md file in folder "day03" which will explain the program (python file).
> is should also include this prompt and the previuos one
