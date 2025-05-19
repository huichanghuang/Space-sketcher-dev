import sys, argparse
import textwrap, importlib
from space_sketcher.__init__ import __version__, __category__
from space_sketcher.tools.text import help_text, sum_help

def safe_import_module(category, pipe):
    try:
        return importlib.import_module("space_sketcher.%s.%s" % (category, pipe))
        # return importlib.import_module("%s.%s" % (category, pipe))
    except ImportError as e:
        print(f"Error: Failed to import 'space_sketcher.{category}.{pipe}'. {e}", file=sys.stderr)
        sys.exit(1)


def safe_getattr(module, attr_name):
    try:
        return getattr(module, attr_name)
    except AttributeError as e:
        print(f"Error: '{module.__name__}' does not have '{attr_name}'. {e}", file=sys.stderr)
        sys.exit(1)

def category_pipe(pipe):
    """
    Retrieves the list of pipes within a given pipe category.
    Imports the '__init__' module for the specified 'pipe', retrieves its '_pipe' attribute, and returns it.

    Args:
        pipe (str): The name of the pipe category.

    Returns:
        list[str]: The list of pipes within the specified category.
    """
    package = importlib.import_module("space_sketcher.%s.__init__"%pipe)
    # package = importlib.import_module("%s.__init__"%pipe)
    pipelist = package._pipe
    return pipelist

def add_pipes_to_subparsers(subparsers, category):
    """
    Adds subparsers for each pipe in a given category.
    Iterates through the pipes in the specified `category`, creates a subparser for each pipe, sets its description
    using `help_sub_text`, and assigns its `func` attribute to the corresponding pipe function.

    Parameters:
    - `subparsers`: `argparse._SubParsersAction` object, the collection of subparsers to which new pipes will be added.
    - `category`: str, the category of pipes to add to the subparsers.
    """
    _pipeList = category_pipe(category)

    for _pipe in _pipeList:
        package = safe_import_module(category, _pipe)
        pipes = safe_getattr(package, _pipe)
        pipes_help = safe_getattr(package, f"helpInfo_{_pipe}")

        parser_step = subparsers.add_parser(
            _pipe,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(f"Help for {category} {_pipe} command")
        )
        pipes_help(parser_step)
        parser_step.set_defaults(func=pipes)

def main():
    """
    Parses command-line arguments and executes the appropriate functionality based on user input.
    Initializes an `argparse.ArgumentParser`, adds global options (e.g., version), and sets up nested subparsers
    for each category and its associated pipes. If no arguments are provided, displays the help message and exits.
    Otherwise, parses the arguments, calls the selected function with the parsed arguments, and exits.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(sum_help())
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=__version__
    )
    subparsers = parser.add_subparsers(dest='parser_step')

    for _category in __category__:
        sub_subparser = subparsers.add_parser(
            _category,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(help_text(_category, 'Space-Sketcher %s' % _category))
        )
        sub_subparser.add_argument(
            '-v',
            '--version',
            action='version',
            version=__version__
        )
        subsub_subparser = sub_subparser.add_subparsers()
        add_pipes_to_subparsers(subsub_subparser, _category)

    if len(sys.argv) == 1 :
        parser.print_help(sys.stderr)
        sys.exit(0)

    elif len(sys.argv) == 2 and sys.argv[1] in __category__: 
        subparsers.choices[sys.argv[1]].print_help(sys.stderr)
        sys.exit(0)
    elif len(sys.argv) == 3 and sys.argv[1] in __category__:
        cat_parser = subparsers.choices[sys.argv[1]]
        subsub_subparser = cat_parser._subparsers._group_actions[0].choices
        if sys.argv[2] in subsub_subparser:
            subsub_subparser[sys.argv[2]].print_help(sys.stderr)
            sys.exit(0)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()