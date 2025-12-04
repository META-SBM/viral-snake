"""Pipeline execution commands"""

import click
from rich.console import Console

console = Console()


@click.group(name='run')
def run_group():
    """🚀 Run pipeline targets"""
    pass


@run_group.command(name='target')
@click.argument('target')
@click.option('--dry-run', is_flag=True, help='Show what would be executed')
def run_target(target, dry_run):
    """Run pipeline for specific target"""
    # TODO: implement
    console.print("[yellow]TODO: Implement pipeline execution[/yellow]")