demo="""<h2>Case study</h2>
<p>In this case study, we are going to perform an analysis of a dataset from corralimorphs on the Palmyra atoll. The study ID of this dataset is 10798, and the full description is available <a href="https://qiita.ucsd.edu/study/description/10798">here</a>. In the original BIOM file, the mapping file is not included. The BIOM file in the <em><strong>massoc </strong></em>data folder has already had the mapping file added.</p>
<h3>Set the default directory and load the data</h3>
<p>In this case, the <em><strong>massoc/data</strong></em> directory was selected. The <strong>Open BIOM files </strong>button will automatically open this directory. Here, we selected the <strong>coral_ID10895_31522.biom </strong>file. As you can see, the BIOM file was imported correctly; the number of samples and taxa is shown in the right text box, as well as the metadata variables.</p>
<p><img src="memory:input.png" alt="Input files for coral case studies" width="500" height="376"/></p>
<h3>Preprocess the file</h3>
<p>Adding a prefix will help us identify the output files, so we filled in <strong>coral</strong>.</p>
<p>Because the dataset is very large, preprocessing steps may help remove some sample heterogeneity and reduce the runtime. If the sample heterogeneity is too large, tools like SPIEC-EASI will not complete network inference successfully.</p>
<p>We removed taxa with a mean count below 2. This will remove singletons and low-abundant taxa. We also rarefy samples to even depth; uneven depths can induce spurious associations. The prevalence filter was set to 40%, because this is a large dataset.</p>
<p>In this case, the variable <strong>taxonomy_string_to_family </strong>contains the host (coral) taxonomy, and could be helpful if we want to construct separate networks for each coral family. Let's select this.</p>
<p>It may also be helpful to construct networks at multiple taxonomic levels. Let's construct networks for the taxonomic units and at the family level.</p>
<p><img src="memory:process.png" alt="Processing settings for coral case studies" width="500" height="378"/></p>
<h3>Select network inference tools and run</h3>
<p>SPIEC-EASI will return sparse, high-precision networks, while CoNet may have higher sensitivity. Both tools are good at identifying hub species. Let's select both and then click <strong>Run network inference. </strong></p>
<p><strong><img src="memory:network.png" alt="Network inference on coral case study" width="500" height="468"/></strong></p>
"""