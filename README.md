# FlowCytometry
## Agnostic Flow cytometry data analyses of peripheral blood mononuclear cells from fibromyalgia patients and controls. 
### Contact Information
* Prepared By: Vivek Verma (vivek.verma@mail.mcgill.ca)

This repository uses a modified version of [VoPo](https://www.nature.com/articles/s41467-020-17569-8>).
The original versions of all the functions is [here](https://github.com/stanleyn/VoPo) 
[opt-SNE](https://github.com/omiq-ai/Multicore-opt-SNE) was implemented using R [wrapper function](https://rdrr.io/github/milescsmith/dim.reduction.wrappers/man/optSNE.html) with minor modifications. 

Following pre-processing steps were performed on FCS files before proceeding with the pipeline:<br />
* Compensated in [FlowJo](https://www.flowjo.com/) v10.6.2<br />
 wsp: S:\FM_FLOW_CYTOMETRY\clean_fcs\1_NK{1..4}.wsp (compensation)<br />
      S:\FM_FLOW_CYTOMETRY\clean_fcs\NK_file_selection.wsp
* Quality control for the flow cytometry data based on compositional data analysis was performed using [flowClean](http://bioconductor.org/packages/release/bioc/html/flowClean.html).
* Good (non-anomalous), alive, singlets were gated and exported using FlowJo:<br />
 wsp: S:\FM_FLOW_CYTOMETRY\clean_fcs\2_NK_gating.wsp
 
 
