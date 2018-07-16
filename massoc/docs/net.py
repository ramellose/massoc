net="""<h2>Network inference</h2>
<p>In this tab, you can specify which tools should be used to run network inference. For some of these tools, you can also adjust the settings. Finally, an overview of the settings is given, as well as an option to copy these to a command line call. </p>
<h3>Select network inference tools to run [-tools]</h3>
<p>Currently, SparCC, CoNet and SPIEC-EASI are supported. When selected, these tools will run network inference using the default settings. </p>
<p>For more details on these tools, please check out the original publications:</p>
<p><a href="http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687">SparCC<br /></a><a href="http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226">SPIEC-EASI<br /></a><a href="https://f1000research.com/articles/5-1519/v1">CoNet</a></p>
<h3>Change settings for SparCC [-spar_boot and -spar_pval]</h3>
<p>SparCC bootstraps correlations, and uses the bootstraps to compute pseudo p-values. Lower p-values can result in sparser networks.</p>
<h3>Change settings for SPIEC-EASI [-spiec]</h3>
<p>SPIEC-EASI can use two different algorithms: the Meinshausen-Buhlmann algorithm and the Graphical Lasso algorithm. For both of these algorithms, the sparsity of the final network is controlled by the number of StARS repetitions.</p>
<p>In the graphical user interface, these settings can be specified directly. For command line, the <strong>-spiec </strong>option allows users to specify an alternative filepath for the SPIEC-EASI RScript (massoc/execs/spieceasi.R). Alternatively, this script can be adjusted to specify the desired algorithm and number of StARS repetitions.</p>
<h3>Change settings for CoNet [-conet]</h3>
<p>As CoNet has a multitude of settings, there is no option to adjust them inside <em><strong>massoc</strong></em>. However, alternative Shell scripts (massoc/execs/CoNet.sh) can be supplied in command line with <strong>-conet</strong>.</p>
<p>The easiest way to change CoNet settings is to generate command line calls from the Cytoscape CoNet app, and to use these to adapt the Shell script. However, care should be taken to preserve the filename format used in the Shell script, as <em><strong>massoc </strong></em>will not be able to read the inferred network otherwise.</p>
<h3>Number of processes [-cores]</h3>
<p>After all settings have been entered, <em><strong>massoc </strong></em>will specify the number of unique networks that need to be inferred. Network inference can be sped up significantly if networks are inferred simultaneously. By default, <em><strong>massoc </strong></em>will distribute the jobs over four cores, but more cores can be added if a particularly large number of networks needs to beinferred.</p>
<h3>Export as command line call</h3>
<p>With all settings entered, <em><strong>massoc </strong></em>can generate a command line call that will perform the same tasks.</p>
<h3>Run network inference</h3>
<p>Network inference will start, and a loading bar keeps track of progress. When finished, one or more output XML files will be written to the previously specified filepath.</p>
"""