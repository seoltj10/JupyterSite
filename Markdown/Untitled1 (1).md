# Interactive volcano plot 

## - Increased secretory genes between wt hepatocyte vs. YAP OFF organoid: microarray data from Yimlamai et al. (2009), Cell


```python
import numpy as np
from bokeh.models import ColumnDataSource, HoverTool, Span, BoxAnnotation
from bokeh.io import push_notebook, show, output_notebook, output_file
from bokeh.plotting import figure, save
import pandas as pd
pd.options.mode.chained_assignment = None
output_notebook()
```



<div class="bk-root">
    <a href="https://bokeh.org" target="_blank" class="bk-logo bk-logo-small bk-logo-notebook"></a>
    <span id="15500">Loading BokehJS ...</span>
</div>





```python
from bokeh.resources import CDN
from bokeh.embed import file_html
```


```python
def volcano(title, total_data, inc_secretome, dec_secretome, normal_color, inc_sec_color, dec_sec_color):
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
    df2 = pd.read_csv(inc_secretome)
    inc_sec = df2["Genesymbol"]                   
    
    df3 = pd.read_csv(dec_secretome)
    dec_sec = df3["Genesymbol"]
    
    not_secreted = df[~(df['Genesymbol'].isin(inc_sec))]
    not_secreted = not_secreted[~(not_secreted['Genesymbol'].isin(dec_sec))]
    not_secreted['Secretion'] = "No"
    
    sig = df[df['adjPVal'] <= 5E-3]
    sig_fc = sig[abs(sig['logFC']) >= 1]
    inc_secreted = sig_fc[sig_fc['Genesymbol'].isin(inc_sec)]
    inc_secreted['Secretion'] = 'Yes'
    dec_secreted = sig_fc[sig_fc['Genesymbol'].isin(dec_sec)]
    dec_secreted['Secretion'] = "Yes"
    
    # Set tooltips; Gene name, logFC value, P value,and secreted or not.                         
    tooltips = [('Gene name', ' @Genesymbol'),('logFC',' @logFC'),('P-value',' @adjPVal'),('Secretion?',' @Secretion')]
    
    # Draw background figure                        
    p = figure(title=title, plot_width=800, plot_height=800, tooltips=tooltips)
    source1 = ColumnDataSource(inc_secreted)
    source2 = ColumnDataSource(dec_secreted)
    source3 = ColumnDataSource(not_secreted)
    
    # Add each circle on figure - overdraw outline upon secreted genes
    p.circle(x="logFC", y="transformed_p", source=source3, size=7, alpha=0.5, legend_label="Not secreted", color=normal_color)
    p.circle(x="logFC", y="transformed_p", source=source1, size=7, alpha=0.8, legend_label='Secreted, Increased', color=inc_sec_color)
    p.circle(x='logFC', y='transformed_p', source=source1, size=7, alpha=1, legend_label='Secreted, Increased', color='black', fill_color=None, name='outlines')
    p.circle(x="logFC", y="transformed_p", source=source2, size=7, alpha=0.8, legend_label='Secreted, Decreased', color=dec_sec_color) 
    p.circle(x='logFC', y='transformed_p', source=source2, size=8, alpha=1, legend_label='Secreted, Decreased', color='black', fill_color=None, name='outlines')

    
    # Add cutoff lines - |logFC| > 1, p value < 5E-4
    p.add_layout(Span(location=1, dimension='height', line_color='black', line_dash='dashed', line_width=1.5))
    p.add_layout(Span(location=-1, dimension='height', line_color='black', line_dash='dashed', line_width=1.5))
    p.add_layout(Span(location=2.30102999566, dimension='width', line_color='black', line_dash='dashed', line_width=1.5))
    
    # Set boxannotation
    left_box = BoxAnnotation(top=2.30102999566, right=-1, fill_alpha=0.3, fill_color='gray')
    right_box = BoxAnnotation(top=2.30102999566, left=1, fill_alpha=0.3, fill_color='gray')
    center_box = BoxAnnotation(left = -1, right=1, fill_alpha=0.3, fill_color='gray')

    p.add_layout(left_box)
    p.add_layout(right_box)
    p.add_layout(center_box)

    # Set titles
    p.title.text_font_size = "13px"
    p.title.align = 'center'
    
    # Set legend properties
    p.legend.border_line_width = 1
    p.legend.border_line_color = "black"
    p.legend.border_line_alpha = 0.5
    p.legend.location = "top_left"
    p.legend.click_policy="hide"
    
    # Set label and background color
    p.xaxis.axis_label = "log2FC"
    p.xaxis.axis_label_text_font_size = "15pt"
    p.yaxis.axis_label = "-log10(P)"
    p.yaxis.axis_label_text_font_size = "15pt"
    p.background_fill_color = "#DFDFE5"
    p.background_fill_alpha = 0.5
    
    # LET IT SHOW ITSELF!                         
    show(p)
    
    # Or save as html file
    output_file('output_plot.html', mode='inline')
    save(p)
```


```python
volcano("Increased secretory genes between wt hepatocyte vs. YAP OFF organoid: microarray data from Yimlamai et al. (2009), Cell",
        "wt hep vs. off org source.csv","wt hep vs. yap off org final value matched -inc.csv","wt hep vs. yap off org final value matched -dec.csv","Green","Red","Blue")
```








<div class="bk-root" id="c5998427-d3a7-423b-8a2e-23952ae96953" data-root-id="14803"></div>




