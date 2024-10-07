volcano_goi = c(
  #TLS related
  # 'TNF', #very low expression
  #'TNFR1', #Receptor for TNF and LTa, removed because produces subscript out of bounds error
  # 'TNFR2', #Receptor for TNF and LTa, , removed because produces subscript out of bounds error
  'TNFSF14', # this is LIGHT
  'LTB',
  'LTA',
  'LTBR', #LT beta receptor, commented out because no difference between responders and non-responders
  'TNFRSF14', #HVEM, commented out because no difference between responders and non-responders
  'TNFRSF17', #BCMA
  # 'TNFRSF8', #CD30, present in activated T and B cells
  'TNFRSF13B', #TACI
  'TNFRSF13C', #BAFF-R  
  'TNFRSF18', #This gene encodes a member of the TNF-receptor superfamily. The encoded receptor has been shown to have increased expression upon T-cell activation, and it is thought to play a key role in dominant immunological self-tolerance maintained by CD25(+)CD4(+) regulatory T cells. Knockout studies in mice also suggest the role of this receptor is in the regulation of CD3-driven T-cell activation and programmed cell death. Three alternatively spliced transcript variants of this gene encoding distinct isoforms have been reported. [provided by RefSeq, Feb 2011] (https://www.genecards.org/cgi-bin/carddisp.pl?gene=TNFRSF18)
  'LAMP3', #DC-lamp
  'GPR183', #EBI-2
  'S1PR4', 
  'CCR7',
  'MADCAM1', 
  # 'CXCL5', 
  'MS4A1', #CD20
  'CD69', 
  'CD79A', #BCR signaling
  'CCL2', #monocyte, DC, and T lymphoctyte trafficking
  'CCL19', 
  'CCL21', 
  'CXCL13',
  'CXCR5', #defining marker for follicular B helper cells (TFH)
  'IL7R', #important for T cell memory formation. IL-7 receptor blockade blunts antigen-specific memory T cell responses and chronic inflammation in primates (https://www.nature.com/articles/s41467-018-06804-y)
  
  'CLEC1B',
  'PDPN',
  'TNFRSF9', #CD137/4-1BB
  'IFNG',
  'GZMK', 
  'CR2', #CD21
  'FCER2', #CD21
  'EOMES',
  'TOX',
  'TCF7',
  'PDCD1',
  'CTLA4',
  'LAG3',
  'HAVCR2', #TIM-3 gene, returns 'subscript out of bounds' when called
  'ENTPD1', #CD39
  #    'SDC1', #CD138, strongly expressed across the board so removed 
  'ICOS',
  'FOXP3',
  'IL2RA', #CD25 / IL2R, marker of Tregs
  'IL2',
  'IL17C', #cytokine reported to stimulate the release of tumor necrosis factor alpha and interleukin 1 beta from a monocytic cell line. The expression of this cytokine was found to be restricted to activated T cells.(https://www.genecards.org/cgi-bin/carddisp.pl?gene=IL17C)
  'IL17REL', #Predicted to enable interleukin-17 receptor activity. Predicted to be involved in cytokine-mediated signaling pathway.(https://www.genecards.org/cgi-bin/carddisp.pl?gene=IL17REL)
  #   'BCL6' #strongly expressed across the board, so removed
  #  'AICDA' #AID gene, uniformly low expression across the board, so removed 
  # 'IFNA1'
  'FCRL1', #Primarily expressed in secondary lymphoid tissues by mature subsets of B-cells (https://www.uniprot.org/uniprot/Q96LA6)
  'FCRL4', #FCRL4 is specifically expressed by memory B cells (MBCs) localized in sub-epithelial regions of lymphoid tissues. Expansion of FCRL4+ B cells has been observed in blood and other tissues in various infectious and autoimmune disorders. (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0179793)
  'FCRL5', #expressed on plasma cells and atypical B cells, induced upon BCR stimulation
  'FCRLA', #FcRLA is uniquely interesting due to its intracellular localization, unusual structural features, and high expression within human germinal center and marginal zone B cells. (https://www.jimmunol.org/content/185/5/2960)
  'IL6',
  'IL22RA1',
  'IL10', #anti-inflammatory cytokine
  # 'HLA-A', 'HLA-B', 'HLA-C',
  # 'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2',
  # 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2'
  # 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5',
  # 'HLA-E', 'HLA-F', 'HLA-G', 'HLA-H', 'HLA-J', 'HLA-L'
  'ADAM28', #important regulator of inflammatory signaling pathways as this protease shed the pro-inflammatory cytokine, pro-TNF-alpha. ADAM28 also interacts with integrins and a P-selectin ligand (PSGL-1) involved in inflammatory cell migration. (https://erj.ersjournals.com/content/52/suppl_62/OA5377)
  'MEF2C', #required for acquisition of effector Treg phenotype and for activation of epigenetic program that suppresses anti-tumor immune responses of conventional T and B cells. Deletion of this gene in Tregs switches off expression of IL10 and ICOS and leads to enhanced antitumor immunity (https://www.frontiersin.org/articles/10.3389/fimmu.2021.703632/full)
  'TRAF3IP3', #essential for T and B cell development. in myeloid cells, regulates the host response to cytosolic viral RNA (https://www.nature.com/articles/s41467-020-16014-0)
  'ELF4', #evidence from human autoinflammatory disease shows that it sustains  the expression of anti-inflammatory genes, such as Il1rn, and limited the upregulation of inflammation amplifiers, including S100A8, Lcn2, Trem1 and neutrophil chemoattractants. Nature Comm 2021 (https://www.nature.com/articles/s41590-021-00984-4)
  'BANK1', #B cell scaffold protein with ankyrin repeats 1, inversely correlated with disease activity in RA mouse model (https://pubmed.ncbi.nlm.nih.gov/29370826/#:~:text=Conclusion%3A%20Decreased%20BANK1%20expression%20promotes,independent%20mechanism%20in%20CIA%20mice); interacts with TRAF6 and MyD88 in innate immune signaling in B cells (https://www.nature.com/articles/s41423-019-0254-9). Per this review excellent article in JEM (https://rupress.org/jem/article/217/10/e20200483/151908/Single-cell-analysis-of-germinal-center-B-cells), BANK1 is a scaffold protein involved in BCR-induced calcium mobilization and in the negative modulation of CD40-mediated AKT activation. BANK1 and RASGRP2 appear to be the earliest specific markers induced in the GC PreM. They found that BANK1 is an early marker of PreM commitment
  'TRAF1', #TNF receptor associated factor. Promotes CD8 memory in LCMV model, deletion leads to decreased memory formation. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3260874/#:~:text=In%20T%20cells%2C%20overexpression%20of,et%20al.%2C%202007). TRAFs regulate the three signals required for the activation, differentiation, and survival of CD4(+) T cells and other T-cell subsets (https://pubmed.ncbi.nlm.nih.gov/26072698/)
  'IL16', #Interleukin 16 (IL-16), formerly known as lymphocyte chemoattractant factor or LCF, is a pro-inflammatory cytokine that is chemotactic for CD4+ T lymphocytes, monocytes, and eosinophils. In addition to inducing chemotaxis, IL-16 can upregulate IL-2 receptor and HLA-DR4 expression, inhibit T cell receptor (TcR)/CD3-dependent activation, and promote repression of HIV-1 transcription. IL-16 is a unique cytokine with no significant sequence homology to other well-characterized cytokines or chemokines.(https://www.bu.edu/interleukin-16/)
  'RIPK2', #T-Cell-Intrinsic Receptor Interacting Protein 2 Regulates Pathogenic T Helper 17 Cell Differentiation (https://www.cell.com/immunity/pdfExtended/S1074-7613(18)30385-6). New therapeutic target in IBD (https://www.frontiersin.org/articles/10.3389/fphar.2021.650403/full)
  'ITGAM', #CD11b, modulates pathogenesis in lupus nephritis (https://www.frontiersin.org/articles/10.3389/fmed.2018.00052/full). α-chain of integrin receptor CD11b/CD18 (also known as αMβ2, Mac-1, and CR3), is highly expressed on the surface of innate immune cells, including macrophages and neutrophils.  CD11b also modulates other signaling pathways in these cells, such as the Toll-like receptor signaling pathways, that mediate generation of type I interferons, a key proinflammatory cytokine and circulating biomarker in SLE and LN patients. 
  'WDFY4', # WDFY4 is required for cross-presentation in response to viral and tumor antigens (https://www.science.org/doi/10.1126/science.aat5030)
  'DOCK10' #Dock10, a novel CZH protein selectively induced by interleukin-4 in human B lymphocytes (https://www.sciencedirect.com/science/article/pii/S0161589008001557). Deletion of Dock10 in B Cells Results in Normal Development but a Mild Deficiency upon In Vivo and In Vitro Stimulations (https://www.frontiersin.org/articles/10.3389/fimmu.2017.00491/full)
  
) 


#### annotations for other significantly DEGs, immuology related
####DDX26B - not well described but does seem to have some relation to CD4 TfH in granulomatosis polyangiitis (https://www.jimmunol.org/content/jimmunol/208/4/807.full.pdf?with-ds=yes)
#IRAK3 - interleukin1 receptor kinase 3, modulates downstream innate immunity signaling (https://www.nature.com/articles/s41598-019-51913-3)
#SLIT2 - axon guidance glycoprotein, modulates inflammatory phenotype in fibrocytes in graves' dz (https://www.jimmunol.org/content/early/2018/05/10/jimmunol.1800259)
#MUC1 - a cancer associated antigen and putative regulatory checkpoint on T cells (https://www.frontiersin.org/articles/10.3389/fimmu.2018.02391/full)
#TNFAIP8L3 - The TIPE (tumor necrosis factor-α-induced protein 8-like) family are newly described regulators of immunity and tumorigenesis consisting of four highly homologous mammalian proteins: TNFAIP8 (tumor necrosis factor-α-induced protein 8), TIPE1 (TNFAIP8-like 1, or TNFAIP8L1), TIPE2 (TNFAIP8L2) and TIPE3 (TNFAIP8L3).(https://www.nature.com/articles/cmi20174)
#MRVI1 a common MRV integration site in BXH2 myeloid leukemias, encodes a protein with homology to a lymphoid-restricted membrane protein Jaw1 (https://www.nature.com/articles/1202419.pdf?origin=ppub). Not very well described...
#CDK17 - not much out there, but has been identified as part of group of gene markers for tolerance after alloBMT https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4361657/
#MMP9 Matrix metalloproteinases (MMP) 3 and 9 as biomarkers of severity in COVID-19 patients (https://www.nature.com/articles/s41598-021-04677-8), MMP9 protects against LPS-induced inflammation in osteoblasts (https://pubmed.ncbi.nlm.nih.gov/31726909/)


#annotattions for other significantly DEGs, not immunology related
#STARD9 - centrosomal protein
#ADAMTSL2  extracellular matrix glycoprotein ADAMTSL2 is increased in heart failure and inhibits TGFβ signalling in cardiac fibroblasts (https://www.nature.com/articles/s41598-021-99032-2)
#DYRK2 - negative regulator of Type I interferon induction (https://pubmed.ncbi.nlm.nih.gov/29155197/)
#PARK15 - regulator of innate immunity http://genesdev.cshlp.org/content/34/5-6/341.full
#ARHGAP15 - negative master regulator of neutrophils (https://pubmed.ncbi.nlm.nih.gov/21551229/), may be FOXP3 target (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5276829/)
#IQSEC2 - associated with genetic disorder a/w autism and intellectual disability
#RAP1GAP2 - GTPase activating protein, function not well characterized
#LUM - collagen binding proteoglycan, involved in ECM
#EPHA3 - receptor tyrosine kinase, possible target in pre-B ALL, published by Soren (https://pubmed.ncbi.nlm.nih.gov/27922598/), thought to be tumor suppressor
#ABI3BP - ECM protein
#ASPN - asporin, ECM protein, involved in cardiac remodelng (https://pubmed.ncbi.nlm.nih.gov/35470068/); described by Ben Ho Park group in restricting mesenchymal stromal cell differentiation (https://aacrjournals.org/cancerres/article/79/14/3636/638218/Asporin-Restricts-Mesenchymal-Stromal-Cell)
#CDH3 - P-cadherin, overexpressed in colorectal CA (https://pubmed.ncbi.nlm.nih.gov/29142905/)
#INMT - methyltransferase, doesn't seem to be well characterized
#SLC34A2 - not a lot, but according to one article may be involved in HCC carcinogenesis, knockdown inhibits HCC (https://www.ingentaconnect.com/content/cog/or/2016/00000024/00000006/art00013;jsessionid=9tgekfm8bioc6.x-ic-live-01)
