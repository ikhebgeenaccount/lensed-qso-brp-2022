import pandas as pd

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


def plots_in_subfigures(gals, plot, caption=True, main_caption_text='', label=None):
    """
    Plots are found in plots/[galaxy]_[plot].pdf
    """
    print('\\begin{figure}\n\centering')
    for i, g in enumerate(gals):
        print('\\begin{subfigure}[t]{.7\\textwidth}\n\centering\n\hspace*{-2.3in}')
        print(f'\includegraphics[width=\linewidth]{{plots/{g}_{plot}.pdf}}')

        if caption:
            print('\caption{' + g + '}')

        if label:
            print(f'\label{{fig:{label}_' + g + '}')
        print('\end{subfigure}%' + ('\n' if i % 2 == 1 else ''))
    if caption:
        print(f'\caption{{{main_caption_text}}}')
    if label:
        print(f'\label{{fig:all_{label}}}')
    print('\end{figure}')


def agnf_output_table(lqsos_df, cols, col_names=None, label=''):
    df_copy = lqsos_df.copy()

    # Table environment: first col is right aligned, rest left
    print('\\begin{table}\n\\caption{}\\label{' + label + '}\\begin{tabular}{|r|' + 'l|' * len(cols) + '}')

    # Create table header with col_names
    header = ''
    for i in range(len(cols)):
        header += (col_names[i] + '&') if col_names else (cols[i] +'&')
    # Exchange last & with newline \\
    header = header[:-1] + '\\\\'

    print(header)
    print('\\hline')  # horizontal line between header and data

    def custom_format(f):
        if f > 1e2:# or f < 1e-3:
            return '{:.2e}'.format(f).replace('e+', 'e')
        else:
            return '{:.2f}'.format(f)

    # Create data table
    lqsos_out = pd.DataFrame()
    for i, c in enumerate(cols):
        has_single_err = has_double_err = False

        if f'{c}_err' in df_copy.columns:
            has_single_err = True
        if f'{c}_pe' in df_copy.columns:
            has_double_err = True

        err_string = ''
        if has_single_err:
            err_string = ' \\pm ' + df_copy[f'{c}_err'].map(custom_format)
        elif has_double_err:
            err_string = '^{+' + df_copy[f'{c}_pe'].map(custom_format) + '}_{-' + df_copy[f'{c}_me'].map(custom_format) + '}'

        is_float = df_copy[c].dtype == float

        if is_float:
            data_col = '$' + df_copy[c].map(custom_format) + err_string + '$'
        else:
            data_col = df_copy[c]

        if i == 0:  # First column
            lqsos_out['out'] = data_col + ' & '
        elif i == len(cols) - 1:  # Last columns
            lqsos_out['out'] += data_col + ' \\\\'
        else:
            lqsos_out['out'] += data_col + ' & '

    for i in range(len(lqsos_out)):
        print(lqsos_out['out'].iloc[i])

    print('\\end{tabular}\n\\end{table}')
