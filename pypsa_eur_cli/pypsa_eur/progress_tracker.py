import re
import sys
from pathlib import Path

from colorama import Fore, Style
from tqdm import tqdm


class ProgressTracker:
    """
    A class to manage progress tracking for PyPSA-Eur workflow execution.

    This class tracks the progress of the workflow execution and provides
    human-readable output for the various rules being executed.
    """

    def __init__(self, rule_action_mapping=None):
        """
        Initialize the progress tracker.

        Args:
            rule_action_mapping: Dictionary mapping rule names to human-readable actions
        """
        self.workflow_progress = None
        self.rule_name = None
        self.printed_statements = set()
        self.rule_action_mapping = rule_action_mapping or {}

        # Regex patterns to match Snakemake output
        self.progress_pattern = re.compile(r"(\d+) of (\d+) steps \((\d+)%\) done")
        # Enhanced patterns to capture more rule types
        self.rule_patterns = [
            re.compile(r"^rule (\w+):"),
            re.compile(r"^localrule (\w+):"),
            re.compile(r"^checkpoint (\w+):"),
            re.compile(r"^localcheckpoint (\w+):"),
        ]

    def format_progress_description(self):
        """Format the progress bar description with consistent color."""
        return f"{Fore.BLUE}Workflow Progress"

    def initialize_progress_bar(self, total_steps):
        """Create progress bar if it doesn't exist yet."""
        if self.workflow_progress is None:
            # Create workflow progress bar
            self.workflow_progress = tqdm(
                total=total_steps,
                desc=self.format_progress_description(),
                unit="step",
                ncols=100,
                dynamic_ncols=True,
                file=sys.stdout,
                position=0,
                bar_format="{desc} {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
            )

    def update_progress(self, completed, total):
        """
        Update progress based on Snakemake output.

        Args:
            completed: Number of completed steps
            total: Total number of steps
        """
        # Initialize progress bar if it doesn't exist yet
        self.initialize_progress_bar(total)

        # Update job count if it changed
        if self.workflow_progress is not None and self.workflow_progress.total != total:
            self.workflow_progress.total = total
            # Re-apply the color formatting to the description
            self.workflow_progress.set_description(self.format_progress_description())

        # Update workflow progress bar
        if self.workflow_progress is not None:
            self.workflow_progress.n = completed
            # Ensure color formatting is maintained after update
            self.workflow_progress.set_description(self.format_progress_description())

    def parse_line(self, line):
        """
        Parse a line from Snakemake output and update progress accordingly.

        Args:
            line: A line of output from Snakemake

        Returns:
            True if the line was handled, False otherwise
        """
        # Look for Snakemake's progress output
        match = self.progress_pattern.search(line)
        if match:
            # Extract progress info
            completed = int(match.group(1))
            total = int(match.group(2))
            self.update_progress(completed, total)
            return True

        # Check for all rule patterns (localrule, rule, checkpoint, etc.)
        for pattern in self.rule_patterns:
            match = pattern.match(line.strip())
            if match:
                rule_name = match.group(1)
                # Immediately print when we detect a rule
                statement_identifier = f"{rule_name}_detected"
                if statement_identifier not in self.printed_statements:
                    action = self.rule_action_mapping.get(
                        rule_name, f'Executing rule "{rule_name}"'
                    )
                    tqdm.write(f"{Fore.CYAN}{action}{Style.RESET_ALL}")
                    self.printed_statements.add(statement_identifier)
                    self.rule_name = rule_name
                return True

        # Also capture rule name from other formats
        if "rule " in line and ":" in line:
            parts = line.split()
            for i, part in enumerate(parts):
                if part in ["rule", "localrule", "checkpoint", "localcheckpoint"] and i + 1 < len(parts):
                    rule_name = parts[i + 1].strip(":")
                    statement_identifier = f"{rule_name}_detected"
                    if statement_identifier not in self.printed_statements:
                        action = self.rule_action_mapping.get(
                            rule_name, f'Executing rule "{rule_name}"'
                        )
                        tqdm.write(f"{Fore.CYAN}{action}{Style.RESET_ALL}")
                        self.printed_statements.add(statement_identifier)
                        self.rule_name = rule_name
                    break

        # Check for lines with wildcards - this provides additional context
        if "wildcards:" in line and self.rule_name:
            # Extract and clean the wildcards string
            wildcard_start = line.find("wildcards: ")
            wildcards_str = line[wildcard_start:].split("resources:")[0].strip()
            wildcards_str = wildcards_str.replace("wildcards:", "").strip()

            # Parse wildcards into a cleaner format
            wildcards_clean = []
            for wc in wildcards_str.split(", "):
                if "=" in wc and not wc.startswith("run="):
                    wildcards_clean.append(wc.strip())

            # Only print wildcard details if we have meaningful wildcards
            if wildcards_clean:
                wildcards_str = ", ".join(wildcards_clean)
                # Create a unique identifier for this specific wildcard combination
                statement_identifier = f"{self.rule_name}_{wildcards_str}"

                if statement_identifier not in self.printed_statements:
                    # Print wildcard details as additional info
                    tqdm.write(f"  {Fore.YELLOW}└─ with {wildcards_str}{Style.RESET_ALL}")
                    self.printed_statements.add(statement_identifier)

            return True

        # Check for job execution messages (for rules without wildcards)
        if "jobid:" in line and self.rule_name:
            # We've already printed the rule when we detected it, so just track it
            statement_identifier = f"{self.rule_name}_jobid"
            self.printed_statements.add(statement_identifier)
            return True

        return False

    def clean_up(self):
        """Clean up progress bars."""
        if self.workflow_progress:
            self.workflow_progress.close()

    def __repr__(self):
        """Return a string representation of the progress tracker state."""
        return f"ProgressTracker()"