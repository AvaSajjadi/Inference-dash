"""
Professor-style import shim:

    from ornor.ModelORNOR import PyModelORNOR

Internally it forwards to the installed nlbayes build.
"""
from nlbayes.ModelORNOR import PyModelORNOR
__all__ = ["PyModelORNOR"]
