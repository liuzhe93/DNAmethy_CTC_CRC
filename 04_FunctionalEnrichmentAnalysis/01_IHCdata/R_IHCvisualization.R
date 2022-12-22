setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/04_FunctionalEnrichmentAnalysis/01_IHCdata")
rm(list = ls())

library("HPAanalyze")

data<-hpaDownload(downloadList = 'histology')
hpaListParam(data)
# $normal_tissue
# $normal_tissue$tissue
# [1] "adipose tissue"    "adrenal gland"     "appendix"          "bone marrow"       "breast"            "bronchus"         
# [7] "caudate"           "cerebellum"        "cerebral cortex"   "cervix"            "colon"             "duodenum"         
# [13] "endometrium 1"     "endometrium 2"     "epididymis"        "esophagus"         "fallopian tube"    "gallbladder"      
# [19] "heart muscle"      "hippocampus"       "kidney"            "liver"             "lung"              "lymph node"       
# [25] "nasopharynx"       "oral mucosa"       "ovary"             "pancreas"          "parathyroid gland" "placenta"         
# [31] "prostate"          "rectum"            "salivary gland"    "seminal vesicle"   "skeletal muscle"   "skin 1"           
# [37] "skin 2"            "small intestine"   "smooth muscle"     "soft tissue 1"     "soft tissue 2"     "spleen"           
# [43] "stomach 1"         "stomach 2"         "testis"            "thyroid gland"     "tonsil"            "urinary bladder"  
# [49] "vagina"            "N/A"               "hypothalamus"      "hair"              "retina"            "lactating breast" 
# [55] "skin"              "thymus"            "cartilage"         "eye"               "pituitary gland"   "dorsal raphe"     
# [61] "choroid plexus"    "substantia nigra"  "sole of foot"     
# 
# $normal_tissue$cell_type
# [1] "adipocytes"                                 "glandular cells"                           
# [3] "lymphoid tissue"                            "hematopoietic cells"                       
# [5] "myoepithelial cells"                        "respiratory epithelial cells"              
# [7] "glial cells"                                "neuronal cells"                            
# [9] "cells in granular layer"                    "cells in molecular layer"                  
# [11] "Purkinje cells"                             "endothelial cells"                         
# [13] "neuropil"                                   "squamous epithelial cells"                 
# [15] "peripheral nerve/ganglion"                  "cells in endometrial stroma"               
# [17] "cardiomyocytes"                             "cells in glomeruli"                        
# [19] "cells in tubules"                           "cholangiocytes"                            
# [21] "hepatocytes"                                "alveolar cells"                            
# [23] "macrophages"                                "germinal center cells"                     
# [25] "non-germinal center cells"                  "ovarian stroma cells"                      
# [27] "exocrine glandular cells"                   "pancreatic endocrine cells"                
# [29] "decidual cells"                             "trophoblastic cells"                       
# [31] "myocytes"                                   "fibroblasts"                               
# [33] "keratinocytes"                              "Langerhans"                                
# [35] "melanocytes"                                "epidermal cells"                           
# [37] "smooth muscle cells"                        "peripheral nerve"                          
# [39] "cells in red pulp"                          "cells in white pulp"                       
# [41] "cells in seminiferous ducts"                "Leydig cells"                              
# [43] "urothelial cells"                           "N/A"                                       
# [45] "follicle cells"                             "chondrocytes"                              
# [47] "elongated or late spermatids"               "pachytene spermatocytes"                   
# [49] "peritubular cells"                          "preleptotene spermatocytes"                
# [51] "round or early spermatids"                  "sertoli cells"                             
# [53] "spermatogonia cells"                        "alveolar cells type I"                     
# [55] "alveolar cells type II"                     "basal cells"                               
# [57] "ciliated cells (cell body)"                 "ciliated cells (cilia axoneme)"            
# [59] "ciliated cells (ciliary rootlets)"          "ciliated cells (tip of cilia)"             
# [61] "goblet cells"                               "cells in basal layer"                      
# [63] "cells in corneal layer"                     "cells in spinous layer"                    
# [65] "extracellular matrix"                       "fibrohistiocytic cells"                    
# [67] "hair follicles"                             "langerhans cells"                          
# [69] "lymphocytes"                                "vascular mural cells"                      
# [71] "Bergmann glia - cytoplasm/membrane"         "Bergmann glia - nucleus"                   
# [73] "GLUC cells - cytoplasm/membrane"            "GLUC cells - nucleus"                      
# [75] "granular cells - cytoplasm/membrane"        "granular cells - nucleus"                  
# [77] "molecular layer - neuropil"                 "molecular layer cells - cytoplasm/membrane"
# [79] "molecular layer cells - nucleus"            "processes in granular layer"               
# [81] "processes in molecular layer"               "processes in white matter"                 
# [83] "Purkinje cells - cytoplasm/membrane"        "Purkinje cells - dendrites"                
# [85] "Purkinje cells - nucleus"                   "synaptic glomeruli - capsule"              
# [87] "synaptic glomeruli - core"                  "white matter cells - cytoplasm/membrane"   
# [89] "white matter cells - nucleus"               "endocrine cells"                           
# [91] "enterocytes"                                "enterocytes - Microvilli"                  
# [93] "mucosal lymphoid cells"                     "glands of Brunner"                         
# [95] "paneth cells"                               "bowman's capsule"                          
# [97] "collecting ducts"                           "distal tubules"                            
# [99] "proximal tubules (cell body)"               "proximal tubules (microvilli)"             
# [101] "arrector pili muscle cells"                 "non-ciliated cells"                        
# [103] "cytotrophoblasts"                           "hofbauer cells"                            
# [105] "syncytiotrophoblasts - cell body"           "syncytiotrophoblasts - microvilli"         
# [107] "sebaceous glands"                           "eccrine glands"                            
# [109] "neuronal projections"                       "synapses"                                  
# [111] "cells in cortex/medulla"                    "cells in cuticle"                          
# [113] "cells in external root sheath"              "cells in internal root sheath"             
# [115] "enterocytes - Gradient"                     "ganglion cells"                            
# [117] "inner nuclear layer"                        "inner plexiform layer"                     
# [119] "limiting membrane"                          "nerve fiber layer"                         
# [121] "outer plexiform layer"                      "photoreceptor cells"                       
# [123] "pigment epithelial cells"                   "lactating glandular cells"                 
# [125] "cells in zona fasciculata"                  "cells in zona glomerulosa"                 
# [127] "cells in zona reticularis"                  "medullary cells"                           
# [129] "sebaceous cells"                            "secretory cells"                           
# [131] "sweat ducts"                                "cortical cells"                            
# [133] "cells in dentate nucleus"                   "corneal epithelial cells"                  
# [135] "hyaloid membrane"                           "lens epithelial cells"                     
# [137] "lens fiber cells"                           "cells in anterior"                         
# [139] "cells in posterior"                         "ductal cells"                              
# [141] "ependymal cells"                           
# 
# 
# $pathology
# $pathology$cancer
# [1] "breast cancer"        "carcinoid"            "cervical cancer"      "colorectal cancer"    "endometrial cancer"  
# [6] "glioma"               "head and neck cancer" "liver cancer"         "lung cancer"          "lymphoma"            
# [11] "melanoma"             "ovarian cancer"       "pancreatic cancer"    "prostate cancer"      "renal cancer"        
# [16] "skin cancer"          "stomach cancer"       "testis cancer"        "thyroid cancer"       "urothelial cancer"

data("hpa_histology_data")
pdf("ColonCancer_PatientPro.pdf")
geneList <- c("TRIP10", "NGFR", "SLC48A1", "SRMS")
cancer_type <- "colorectal cancer"
hpaVisPatho(data=hpa_histology_data,
            targetGene=geneList,
            targetCancer = cancer_type)
dev.off()


pdf("Colon_TissueORCell.pdf")
geneList <- c("TRIP10", "NGFR", "SLC48A1", "SRM")
tissueList <- "colon"
hpaVisTissue(data=hpa_histology_data,
            targetGene=geneList,
            targetTissue = tissueList)
dev.off()

