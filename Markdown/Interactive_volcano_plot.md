
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
import numpy as np
from bokeh.models import ColumnDataSource, HoverTool, Span
from bokeh.io import push_notebook, show, output_notebook, output_file
from bokeh.plotting import figure, save
import pandas as pd
pd.options.mode.chained_assignment = None
output_notebook() # you can omit this line if you're not using jupyter notebook
```


```python
# These lines are used for generating html file
from bokeh.resources import CDN
from bokeh.embed import file_html
```

 I made a funtion which uses 5 parameters - 

1) Title, 2) Name of total data CSV file to draw background dots, 3) Name of secretome CSV file to pick up secreted genes within whole dataset, 4) Name of color for backgorund genes, and 5) Name of color for secreted genes. those are all string type.

**Make sure upload your CSV files within a same folder with notebook ipynb file!!**


```python
# This block is used for upload CSV file from your local drive to Colab environment - For first CSV
from google.colab import files
uploaded = files.upload()
```



     <input type="file" id="files-d2c700c3-a60b-46bc-a5de-cd9aeb1faf42" name="files[]" multiple disabled />
     <output id="result-d2c700c3-a60b-46bc-a5de-cd9aeb1faf42">
      Upload widget is only available when the cell has been executed in the
      current browser session. Please rerun this cell to enable.
      </output>
      <script src="/nbextensions/google.colab/files.js"></script> 


    Saving 1wk_final value matched.csv to 1wk_final value matched.csv
    


```python
# This block is used for upload CSV file from your local drive to Colab environment - For second CSV
from google.colab import files
uploaded = files.upload()
```



     <input type="file" id="files-b57b86ff-a754-4914-885c-69a03b986a1b" name="files[]" multiple disabled />
     <output id="result-b57b86ff-a754-4914-885c-69a03b986a1b">
      Upload widget is only available when the cell has been executed in the
      current browser session. Please rerun this cell to enable.
      </output>
      <script src="/nbextensions/google.colab/files.js"></script> 


    Saving volcano_1wk_for data.csv to volcano_1wk_for data.csv
    


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
    p.circle(x="logFC", y="transformed_p", source=source2, size=4, color=normal_color, alpha=0.5, legend_label="Not secreted gene")
    p.circle(x="logFC", y="transformed_p", source=source1, size=4, color=secreted_color, alpha=0.5, legend_label="Secreted gene")
    p.circle(x='logFC', y='transformed_p', size=5,alpha=1,source=source1, color='black', fill_color=None, name='outlines', legend_label="Secreted gene")

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










  <div class="bk-root" id="2368fd32-fdc5-493d-b3bd-ccfac9059a36" data-root-id="5488"></div>




