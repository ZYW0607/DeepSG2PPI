# DeepSG2PPI

DeepSG2PPI: A Protein-Protein Interaction Prediction Method Based on Deep Learning Binding Sequence and Functional Features.

A three-input neural network integrating 1D-CNN and 2D-CNN combined with coordinate attention(CA) mechanism is constructed for PPI prediction.


# Requirement

Python==3.6.13   
keras==2.4.3   
numpy==1.19.5   
tensorflow-gpu==2.4.0

GPU acceleration uses CUDA11.2 version, and the corresponding cuDNN version is 8.1.


# Data needed to train the model

Download URL of protein sequence and PPI network data:https://cn.string-db.org/cgi/download.

The relationship file download URL between GO terms:http://geneontology.org/docs/download-ontology/.

Data download link between protein and GO term:https://www.uniprot.org/uniprotkb.


# Train with your own data

If you want to use this network to train your own model, you need to prepare a file containing protein interactions (each line in the file is an interacting protein pair, and the two proteins are separated by a space), a file containing protein ids (file One protein id per line in the file), a file containing protein sequences (one protein sequence per line in the file), where the lines in the protein id file should be relative to the lines in the protein sequence file (that is, the protein sequence file The sequence in line ith corresponds to the protein in line ith of the protein ID file).First use get_npy.py to rewrite the data into .npy format files, and then use Sequence_coding.py to encode protein sequence and local context information, and global statistical information of protein sequence.

In addition, the relationship data between proteins and GO functions needs to be prepared, and the specific relationship preprocessing method can be obtained from the DeepSG2PPI paper.Graph embedding vector representations of proteins were obtained from their graphs against GO function using a node2vec model.The code of the node2vec method has been put into the node2vec_src folder of this repository. For the detailed usage method, please refer to: https://github.com/aditya-grover/node2vec. Then use get_GO_features.py to obtain the GO functional features of the protein through the obtained graph embedding vector. The detailed encoding process can be obtained through the DeepSG2PPI paper.

After the data is encoded, three .npy files will be obtained (respectively, protein sequence and local context information, global statistical information of protein sequence, and graph embedding vector features), which are used as the input of the network.Use the model.py file in the DeepSG2PPI_model folder to train the PPI prediction model.In model.py, the positive and negative samples will be divided and given corresponding labels. The specific division ratio can be adjusted by the test_rate parameter in the file.Of course, you can also add performance indicators you feel credible to the file (if you add the corresponding evaluation indicators yourself, when you use the trained model for prediction, you need to make corresponding predictions in the prediction file PPI_forecast.py Revise).


# Predict your data based on the trained model

The trained model Trained_model.h5 is uploaded in this repository.If you want to use the model to predict your own test data, you must prepare the test data as a .txt (or .csv) file with one proteinName per line.

In the Predict folder, the PPI_forecast.py file is then run to predict whether all possible protein pair combinations in the file have interactions. Of course, you can also predict specific protein pairs with PPI_forecast_p-p.py. The prediction results will output a protein pair file PP_interaction.txt with interactions (where each line represents a protein pair with interactions, and the two proteins are separated by tabs), and a protein pair file PP_noninteraction.txt without interactions ( Each row represents a protein pair that has no interaction, and the two proteins are separated by the same tab).

# Cite

Please cite our paper if you use this code in your own work:

Fan Zhang, Yawei Zhang, Xiaoke Zhu, Xiaopan Chen, Fuhao Lu, and Xinhong Zhang. DeepSG2PPI: A protein-protein interaction prediction method based on deep learning. IEEE/ACM Transactions on Computational Biology and Bioinformatics, pages 1–14, 2023. https://doi.org/10.1109/TCBB.2023.3268661

@article{zhang2023DeepSG2PPI,

  author = {Zhang, Fan and Zhang, Yawei and Zhu, Xiaoke and Chen, Xiaopan and Lu, Fuhao and Zhang, Xinhong},
  
  year ={2023},
  
  pages ={1--14},
  
  title = {DeepSG2PPI: A Protein-Protein Interaction Prediction Method Based on Deep Learning},
  
  journal ={IEEE/ACM Transactions on Computational Biology and Bioinformatics},
  
  doi = {https://doi.org/10.1109/TCBB.2023.3268661}
  
}

# Contact

If you need any assistance, please feel free to contact us:zhangyawei@henu.edu.cn
