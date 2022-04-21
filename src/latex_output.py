NOTHING_FORMATTER = lambda s: str(s)
FLUX_FORMATTER = lambda s: '{:.3e}'.format(s) if s != 0. else '-'

COL_TITLES = {
    'wavelength': '$\lambda (\\unit{\\angstrom})$',
    'filter': 'Filter',
    'source': 'Source',
    'upper_limit': 'Upper limit?',
    'telescope': 'Telescope'
}

def _get_header(col):
    if 'flux' not in col:
        return COL_TITLES[col]
    else:
        comp = col.split('_')[1]

        # Special case for flux_err
        if comp == 'err':
            comp = 'total'

        if 'err' in col:
            return f'$\sigma_{{F_\mathrm{{{comp}}}}} (\\unit{{\milli\jansky}})$'
        else:
            return f'$F_\mathrm{{{comp}}} (\\unit{{\milli\jansky}})$'


def _escape_underscore_str_cols(df):
    # Find columns that have str datatype
    str_cols = df.dtypes[df.dtypes == 'object'].index.tolist()
    # Replace _ with \_
    for c in str_cols:
        df[c] = df[c].str.replace('_', '\_')

    return df


def dataframe_to_latex_table(dfr, outline=None, usecols=None, label=None, caption=None, header=True, formatters=None, **kwargs):
    df = dfr.copy()

    _escape_underscore_str_cols(df)

    print(df.to_latex(columns=usecols, index=False, label=label, caption=caption, header=header, escape=False, formatters=formatters, **kwargs))


def sed_to_latex_table(lqso):
    print('\\begin{landscape}')
    print(f'\\subsection{{{lqso.name}}}')
    df = lqso.filter_sed()

    formatters = []
    usecols = []
    header = []

    for col in df.columns:
        if 'flux' in col and sum(df[col]) > 0 and 'demag' not in col:
            formatters.append(FLUX_FORMATTER)
            usecols.append(col)
            header.append(_get_header(col))
        elif col in COL_TITLES:
            formatters.append(NOTHING_FORMATTER)
            usecols.append(col)
            header.append(_get_header(col))

    dataframe_to_latex_table(df, usecols=usecols, label=f'sed_table:{lqso.name}', formatters=formatters, header=header, position='h!')
    print('\\end{landscape}')


def plots_in_subfigures(gals, plot, caption):
    """
    Plots are found in plots/[galaxy]_[plot].pdf
    """
    print('\\begin{figure}\n\centering')
    for i, g in enumerate(gals):
        print('\\begin{subfigure}{.5\\textwidth}\n\centering')
        print(f'\includegraphics[width=\linewidth]{{plots/{g}_{plot}.pdf}}')
        # print('\caption{' + g + '}\n\label{fig:sed_' + g + '}')
        print('\end{subfigure}%' + ('\n' if i % 2 == 1 else ''))
    print(f'\caption{{{caption}}}\n\label{{fig:all_seds}}\n\end{{figure}}')


def all_seds_plot(gals):
    plots_in_subfigures(gals, 'SED_total', 'Spectral energy distributions (SEDs) of all galaxies.')



