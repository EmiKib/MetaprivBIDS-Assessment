# MetaprivBIDS-Assessment
Statistical assessment of  MetaprivBIDS https://github.com/CPernet/metaprivBIDS. 

Assessed datatset: 

1. AOMIC-ID1000-ds003097<br>
link: https://openneuro.org/datasets/ds003097/versions/1.2.1

2. EEG-Study-Alcohol-Imagery-Task-ds004515<br>
link: https://openneuro.org/datasets/ds004515/versions/1.0.0

3. Dallas-Lifespan-Brain-Study-ds004856<br>
link: https://openneuro.org/datasets/ds004856/versions/1.2.0

4. NIMH-Healthy-Research-Volunteer-Data-ds004215<br>
link: https://openneuro.org/datasets/ds004215/versions/2.0.1

5. The-Midnight-Scan-Club-dataset-ds000224<br>
link: https://openneuro.org/datasets/ds000224/versions/1.0.4

6. Stressful experiences are associated with reduced neural responses to naturalistic emotional and social content in children ds004228<br>
link: https://openneuro.org/datasets/ds004228/versions/1.0.1/file-display/phenotype:questionnaires.json

7. Multivariate Assessment of Inhibitory Control in Youth: Links with Psychopathology and Brain Function Dataset ds004935<br>
link: https://openneuro.org/datasets/ds004935/versions/1.0.2

### Row level correlation:

| Dataset                 | Pearson Correlation (SUDA vs PIF) | Spearman Correlation (SUDA vs PIF) | Kendall's Tau (SUDA vs PIF)  |
|-------------------------|-----------------------------------|------------------------------------|------------------------------|
| AOMIC ds003097          | 0.74                              | 0.76                               | 0.61                         |
| DALLAS ds004856         | 0.57                              | 0.33                               | 0.40                         |
| EGG ds004515            | 0.52                              | 0.68                               | 0.56                         |
| MIDNIGHT ds000224       | 0.74                              | 0.69                               | 0.63                         |
| NIMH ds004215           | 0.50                              | 0.83                               | 0.63                         |
| STRESSFUL ds004228      | 0.47                              | 0.82                               | 0.63                         |
| MULTIVARIATE ds004935   | na                                | na                                 | na                           |
| 8                       | na                                | na                                 | na                           |
| 9                       | na                                | na                                 | na                           |
| 10                      | na                                | na                                 | na                           |


### Field level correlation:



#### Dataset: AOMIC (Field Level)

| Correlation Type           | Variables             | Correlation Value |
|----------------------------|-----------------------|--------------------|
|**Pearson Correlation**    | SUDA & PIF            | 0.77              |
|                            | K-combined & PIF      | 0.37              |
|                            | SUDA & K-combined     | -0.09             |
|**Spearman Rank Correlation** | PIF & SUDA       | 0.82              |
|                            | PIF & K-combined               | 0.11              |
|                            | SUDA & K-combined              | -0.01             |




###  Sum of combined risk for all combination of fields at a minimum of 3 variables, correlation:
