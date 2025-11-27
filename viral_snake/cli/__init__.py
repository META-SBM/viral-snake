"""Command-line interface for viral-snake"""

import click
from rich.console import Console

from . import info, validate, create, run

console = Console()


@click.group()
@click.version_option()
@click.pass_context
def main(ctx):
    """🐍 Viral-Snake: Viral metagenomics pipeline"""
    ctx.ensure_object(dict)


# Register command groups
main.add_command(info.info_group)
main.add_command(validate.validate_group)
main.add_command(create.create_group)
main.add_command(run.run_group)


if __name__ == '__main__':
    main()