
# Draw an interactive volcano plot using bokeh library


#### Even there are so many excellent libraries to visualize data, but bokeh provides easy way to draw **"interactive"** plot!

written by Taejun Seol, 2020-04-18

Here are packages you need to install and import... if you don't have these. try to install them using 

    pip install *package name*

or

    conda install *package name*

If you are using jupyter notebook, you only need to install bokeh. So, just type

    !pip install bokeh

(**! in front of pip is essential**)


```python
!pip install bokeh
```

    Requirement already satisfied: bokeh in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (1.2.0)
    Requirement already satisfied: tornado>=4.3 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from bokeh) (6.0.3)
    Requirement already satisfied: numpy>=1.7.1 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from bokeh) (1.16.4)
    Requirement already satisfied: six>=1.5.2 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from bokeh) (1.12.0)
    Requirement already satisfied: Jinja2>=2.7 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from bokeh) (2.10.1)
    Requirement already satisfied: packaging>=16.8 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from bokeh) (19.0)
    Requirement already satisfied: python-dateutil>=2.1 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from bokeh) (2.8.0)
    Requirement already satisfied: PyYAML>=3.10 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from bokeh) (5.1.1)
    Requirement already satisfied: pillow>=4.0 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from bokeh) (6.1.0)
    Requirement already satisfied: MarkupSafe>=0.23 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from Jinja2>=2.7->bokeh) (1.1.1)
    Requirement already satisfied: pyparsing>=2.0.2 in c:\users\seoltj10\appdata\local\continuum\anaconda3\lib\site-packages (from packaging>=16.8->bokeh) (2.4.0)
    


```python
import numpy as np
from bokeh.models import ColumnDataSource, HoverTool, Span
from bokeh.io import push_notebook, show, output_notebook, output_file
from bokeh.plotting import figure, save
import pandas as pd
pd.options.mode.chained_assignment = None
output_notebook() # you can omit this line if you're not using jupyter notebook
```



    <div class="bk-root">
        <a href="https://bokeh.pydata.org" target="_blank" class="bk-logo bk-logo-small bk-logo-notebook"></a>
        <span id="1001">Loading BokehJS ...</span>
    </div>





```python
# thses lines are used for generating html file
from bokeh.resources import CDN
from bokeh.embed import file_html
```

 I made a funtion which uses 5 parameters - 

1) Title, 2) Name of total data CSV file to draw background dots, 3) Name of secretome CSV file to pick up secreted genes within whole dataset, 4) Name of color for backgorund genes, and 5) Name of color for secreted genes. those are all string type.

**Make sure upload your CSV files within a same folder with notebook ipynb file!!**


```python
# Here is the function to process data and draw plot
def volcano(title, total_data, secretome, normal_color, secreted_color):
    # Start to process whole differential gene expression data to draw plot
    df1 = pd.read_csv(total_data)
    df = df1.copy()
    
    # Remove dot in column name(dots can not be recognized by tooltip functions)
    if "Gene.symbol" in df.columns:
        df.rename(columns = {'Gene.symbol' : 'Genesymbol'}, inplace = True)
    if "adj.P.Val" in df.columns:
        df.rename(columns = {'adj.P.Val' : 'adjPVal'}, inplace = True)
    
    # Transfrom p values into -log10(p)                        
    transformed_p = -df["adjPVal"].apply(np.log10)
    df['transformed_p'] = transformed_p
    
    # Start to process secretome lists - discriminate secretome with others                   
    df2 = pd.read_csv(secretome)
    secretome=df2["Genesymbol"]                   
    not_secreted = (df[~df['Genesymbol'].isin(secretome)])
    not_secreted['Secretion'] = "No"
    sig = df[df['adjPVal'] <= 5E-4]
    sig_fc = sig[sig['logFC'] >= 1]
    is_secreted = sig_fc[sig_fc['Genesymbol'].isin(secretome)]
    is_secreted['Secretion'] ='Yes'
    
    # Set tooltips; Gene name, logFC value, P value,and secreted or not.                         
    tooltips = [('Gene name', '@Genesymbol'),('logFC','@logFC'),('P-value','@adjPVal'),('Secretion?','@Secretion')]
    
    # Draw background figure                        
    p = figure(title=title, plot_width=900, plot_height=900, tooltips=tooltips)
    source1 = ColumnDataSource(is_secreted)
    source2 = ColumnDataSource(not_secreted)
    
    # Add each circle on figure - overdraw outline upon secreted genes
    #NOTE: in jupyter notebook, legend_labels can not be reconginzed as correc attribute, use legend instead.
    p.circle(x="logFC", y="transformed_p", source=source2, size=5, legend="Not secreted", color=normal_color, alpha=0.5)
    p.circle(x="logFC", y="transformed_p", source=source1, size=5, legend='Secreted', color=secreted_color, alpha=0.8)
    p.circle(x='logFC', y='transformed_p', size=6,alpha=1,source=source1, legend='Secreted', color='black', fill_color=None, name='outlines')

    # Add cutoff lines - |logFC| > 1, p value < 5E-4
    p.add_layout(Span(location=1, dimension='height', line_color='black', line_dash='dashed', line_width=1.5))
    p.add_layout(Span(location=-1, dimension='height', line_color='black', line_dash='dashed', line_width=1.5))
    p.add_layout(Span(location=3.30102999566, dimension='width', line_color='black', line_dash='dashed', line_width=1.5))

    # Set titles
    p.title.text_font_size = "18px"
    p.title.align = 'center'

    
    # Set label and background color
    p.xaxis.axis_label = "log2FC"
    p.yaxis.axis_label = "-log10(Pval)"
    p.background_fill_color = "#DFDFE5"
    p.background_fill_alpha = 0.5
    
    # LET IT SHOW ITSELF!                         
    show(p)
    
    # Or save as html file
    output_file('output_plot.html', mode='inline')
    save(p)
    
    
volcano("Increased secretory genes in microarray data from Yimlamai et al. (2009), Cell",
        "volcano_1wk_for data.csv","1wk_final value matched.csv","Green","Red")
```








  <div class="bk-root" id="93a9ba83-e806-4628-a281-0638fd6ec723" data-root-id="2037"></div>






```python

```
