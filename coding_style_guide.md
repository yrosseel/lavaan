---
editor_options: 
  markdown: 
    wrap: 72
---

# Welcome to lavaan coding style guide <!-- omit in toc -->

The guidelines below are primarily intended for new code. Existing code has
already been partially modified, but modification isn't always possible or
desirable because dependent packages and frequently used functions must remain
functional, or because certain names are written in a different style
throughout the literature.

## Names for functions 

For exported functions the names of functions other then `lavaan`, `sem`,
`sam`, `cfa`, `efa` and `growth` should start with 'lav' and use Camel Case,
e.g.  `lavPredict`.

For internal functions the names of functions should start with 'lav\_' and be
in Snake Case, e.g. `lav_fit_cfi_lavobject`.

## Names for function arguments and variables.

The names for function arguments and variables should be lowercase and dotted,
e.g. `categorical.flag`.

## Warnings, errors and notes.

To have a standard way of showing messages one should never use the R functions
`message`, `warning` or `stop` directly, but call the (lavaan internal)
functions `lav_msg_note`, `lav_msg_warn` or `lav_msg_stop` instead.

## Formatting messages

Message should always be formatted with the R functions `gettext`, `gettextf`
or `ngettext` to make them suitable for translation. 

- All messages should begin with a capital letter and end with ".".
- All function names, argument names, and argument values (except numeric 
  literals) are passed as parameters to `gettextf`.
- Lists of values are formatted with internal function `lav_msg_view`. 
- When referring to arguments in a function, use the `=` symbol after the
  argument name to make clear it is an argument; for example, when the message
  concerns the `data` argument in the function `sem()`, the message could
  start with `data= argument is not a data.frame, ...`
