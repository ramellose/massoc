input="""<h2>Input files</h2>
<p>In this tab, you can supply files to <em>massoc</em>. Moreover, you can save or clear the current settings, and load saved settings from a .xyz settings file. If you supplied the correct inputs, the <strong>dialog box</strong> on the right will summarize some properties of your files. </p>
<h3>Set default directory [-fp]</h3>
<p>This directory contains the filepath to a location where you would like output files to be stored.</p>
<h3>Open BIOM files [-biom]</h3>
<p>Specify the complete filepaths of your BIOM files. You can select multiple BIOM files at once, although you can only specify <em>massoc</em>'s settings once; hence, if these settings are incompatible with some of your BIOM files, this may cause <em>massoc </em>to crash. To make sure that <em>massoc </em>works, only run one BIOM file at a time or make sure that your BIOM files have similar metadata.</p>
<h3>Open count tables [-otu]</h3>
<p>Specify the complete filepaths of your tab-delimited count tables. You can select multiple files at once, but make sure these are selected in the same order as your taxonomy files. The count tables should be in a format that can be accept by the BIOM-format parser. This looks like the table below, with a # as the first character, sample identifiers as columns and taxon identifiers as rows:</p>
<p>&nbsp;</p>
<table style="height: 100px;" border="2" width="444">
<tbody>
<tr>
<td style="width: 104px;">#OTU</td>
<td style="width: 104px;">sample_1</td>
<td style="width: 104px;">sample_2</td>
<td style="width: 104px;">sample_3</td>
</tr>
<tr>
<td style="width: 104px;">otu_1</td>
<td style="width: 104px;">7</td>
<td style="width: 104px;">334</td>
<td style="width: 104px;">5</td>
</tr>
<tr>
<td style="width: 104px;">otu_2</td>
<td style="width: 104px;">0</td>
<td style="width: 104px;">189</td>
<td style="width: 104px;">27</td>
</tr>
</tbody>
</table>
<h3>Open taxonomy tables [-tax]</h3>
<p>Specify the complete filepaths of your tab-delimited taxonomy tables. The format for these tables is very similar to the count tables. Again, the # is important for the parser to recognize the column headers. </p>
<p>&nbsp;</p>
<table style="height: 100px; width: 750px;" border="2">
<tbody>
<tr>
<td style="width: 104px;">#OTU</td>
<td style="width: 104px;">Kingdom</td>
<td style="width: 104px;">Phylum</td>
<td style="width: 104px;">Class</td>
<td style="width: 104px;">Order</td>
<td style="width: 104px;">Family</td>
<td style="width: 104px;">Genus</td>
<td style="width: 104px;">&nbsp;Species</td>
</tr>
<tr>
<td style="width: 104px;">otu_1</td>
<td style="width: 104px;">&nbsp;Bacteria</td>
<td style="width: 104px;">&nbsp;Proteobacteria</td>
<td style="width: 104px;">&nbsp;Alphaproteobacteria</td>
<td style="width: 104px;">&nbsp;Rhizobiales</td>
<td style="width: 104px;">&nbsp;...</td>
<td style="width: 104px;">&nbsp;...</td>
<td style="width: 104px;">&nbsp;...</td>
</tr>
<tr>
<td style="width: 104px;">otu_2</td>
<td style="width: 104px;">&nbsp;Bacteria</td>
<td style="width: 104px;">&nbsp;Acidobacteria</td>
<td style="width: 104px;">&nbsp;Acidobacteria</td>
<td style="width: 104px;">&nbsp;Acidobacteriales</td>
<td style="width: 104px;">&nbsp;...</td>
<td style="width: 104px;">&nbsp;...</td>
<td style="width: 104px;">&nbsp;...</td>
</tr>
</tbody>
</table>
<h3>Open metadata files [-tax]</h3>
<p>Specify the complete filepaths of your tab-delimited metadata tables. These should have the sample identifiers as row names.&nbsp;</p>
<p>&nbsp;</p>
<table style="height: 100px;" border="2" width="444">
<tbody>
<tr>
<td style="width: 104px;">#ID</td>
<td style="width: 104px;">pH</td>
<td style="width: 104px;">altitude</td>
<td style="width: 104px;">kelvin</td>
</tr>
<tr>
<td style="width: 104px;">otu_1</td>
<td style="width: 104px;">6.23</td>
<td style="width: 104px;">4538</td>
<td style="width: 104px;">281.2</td>
</tr>
<tr>
<td style="width: 104px;">otu_2</td>
<td style="width: 104px;">5.49</td>
<td style="width: 104px;">3780</td>
<td style="width: 104px;">286.8</td>
</tr>
</tbody>
</table>
<h3>Save, load or clear settings</h3>
<p>When you want to save the current settings, first specify the settings in the other tabs and then return here. The S<strong>ave settings </strong>button saves the settings to a .xyz file. This file can also be opened in Notepad for manual editing. Each row corresponds to a particular <strong><em>massoc </em></strong>setting. The <strong>Load settings </strong>button will set <em><strong>massoc</strong></em>'s settings to the settings in the selected file and display these settings in the graphical user interface. The <strong>Clear settings </strong>button removes all current settings from <em><strong>massoc</strong></em>.</p>
<p>&nbsp;</p>
"""