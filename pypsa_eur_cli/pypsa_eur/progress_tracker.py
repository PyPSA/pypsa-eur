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
        
        # Regex pattern to match Snakemake progress output
        self.progress_pattern = re.compile(r"(\d+) of (\d+) steps \((\d+)%\) done")
    
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
        
        # Look for the 'localrule' to capture the rule name
        if "localrule" in line:
            self.rule_name = line.split()[1].strip(":")  # Remove the colon if present
            return True
        
        # Check for lines with 'localcheckpoint' (special case for checkpoints)
        elif "localcheckpoint" in line:
            checkpoint_name = line.split()[1].strip(":")
            # Get the action for the checkpoint
            action = self.rule_action_mapping.get(
                checkpoint_name, f'Running checkpoint "{checkpoint_name}"'
            )
            
            # Print checkpoint action
            if checkpoint_name not in self.printed_statements:
                tqdm.write(f"{action}")
            self.printed_statements.add(checkpoint_name)
            return True
        
        # Check for lines with wildcards
        elif "wildcards:" in line:
            # Extract and clean the wildcards string
            wildcard_start = line.find("wildcards: ")
            wildcards_str = line[wildcard_start:].split("resources:")[0].strip()
            wildcards_str = wildcards_str.replace("wildcards:", "").strip()
            
            # Create a unique identifier for the rule and wildcards
            statement_identifier = f"{self.rule_name} ({wildcards_str})"
            
            # Get the action from the rule action mapping, default to "Running rule"
            action = self.rule_action_mapping.get(
                self.rule_name, f'Running rule "{self.rule_name}" for'
            )
            
            # Print rule with action and wildcards, only if not already printed
            if statement_identifier not in self.printed_statements:
                tqdm.write(f"{action} {wildcards_str}")
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