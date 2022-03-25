

def dataframe_to_latex_table(dfr, outline=None, usecols=None, label=None, caption=None, header=True):
    df = dfr.copy()

    if usecols is None:
        usecols = df.columns

    # Find columns that have str datatype
    str_cols = df.dtypes[df.dtypes == 'object'].index.tolist()
    # Replace _ with \_
    for c in str_cols:
        if c in usecols:
            df[c] = df[c].str.replace('_', '\_')

    print(df.to_latex(columns=usecols, index=False, label=label, caption=caption, header=header, escape=False))
