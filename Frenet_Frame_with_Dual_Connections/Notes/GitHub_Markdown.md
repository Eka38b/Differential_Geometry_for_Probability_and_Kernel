# GitHub README math conventions

All README and note files in this archive follow these conventions.

1. Use GitHub inline math with dollar delimiters and display equations with double-dollar delimiters on separate lines.
2. Avoid the operator-name macro rejected by the target renderer. Prefer built-in commands such as `\det`, `\min`, and `\sum`, or ordinary roman text with `\mathrm{...}`.
3. Put every subscript in braces, including single indices. Write `$N_{i}^{*}$`, `$D_{i,i+1}$`, `$S_{\perp}$`, and `$e_{1}$`.
4. For a named indexed object, put its whole name or index inside one brace group, or use a superscript for a label. Do not append a second subscript to an already subscripted expression without grouping it.
5. Keep labels and long explanations outside equations where possible. Avoid custom macros and renderer-specific extensions.

The archive check validates local Markdown links, brace structure, absence of the rejected macro, and bracing of subscripts. Optional LaTeX compilation is attempted only when a complete local installation is available. These are source checks; they do not claim that a live GitHub page was opened and inspected.
