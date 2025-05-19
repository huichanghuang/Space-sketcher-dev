import typer
from typing_extensions import Annotated
from types import SimpleNamespace
from enum import Enum
from importlib import metadata

_version = metadata.version("space_sketcher")


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False

def version_callback(value: bool):
    if value:
        print(f"{_version}")
        raise typer.Exit()

app = typer.Typer(context_settings={"help_option_names": ["-h", "--help"]},help=f"Spatial Transcrptomic Analysis ({_version})",)

@app.callback()
def version(
    version: bool = typer.Option(None, '-v',"--version", callback=version_callback, is_eager=True,help="show version."),
):
    return

@app.command(help="Spatial Transcrptomic Analysis")
def count():
    from space_sketcher.rna.run import run
    run()

def main():
    app()

if __name__ == "__main__":
    main()