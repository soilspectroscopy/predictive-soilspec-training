# Introduction

**Predictive soil spectroscopy (PSS)** is a technique that explores the interaction of electromagnetic radiation with soils to estimate several properties via machine learning or chemometrics models. It is a powerful tool for characterizing, measuring, and monitoring soils, as it can provide rapid, cost-effective, and environmentally benign measurements for a number of soil properties from a single scan. However, it has some aspects that require careful consideration during project development, such as:

- **Fitness-for-purpose**: the choice of methodology and spectral range of interest should fit a specific objective and be suitable for the study or application in consideration.  
- **Good predictions flow from good data**: the quality and representativeness of the reference [soil spectral library (SSL)](https://www.woodwellclimate.org/open-soil-spectral-library/), including analytical reference data, are crucial for the accuracy and reliability of the predictions.  
- **Good practices for model building**: the selection and evaluation of the multivariate calibration or machine learning algorithm, and the validation and interpretation of the PSS models and predictions are essential to ensure the robustness and applicability of results.  

In predictive soil spectroscopy, we rely on a reference SSL — a collection of soil spectra and corresponding soil properties measured by conventional methods. The SSL serves as the training dataset for building and testing predictive models. Its performance depends largely on how well it represents the variability and diversity of the soil samples we want to predict. We therefore have two options: create a new SSL tailored to a specific goal, or use an existing public SSL that is sufficiently representative for our analysis. Examples of such libraries are the [**USDA NRCS Kellogg Soil Survey Laboratory Spectral Library**](https://ncsslabdatamart.sc.egov.usda.gov/) and the [**Open Soil Spectral Library**](https://soilspectroscopy.github.io/ossl-manual/).

## Fitness-for-purpose

Soil spectroscopy is a fit-for-purpose technology, which means that the predictive solution must be suitable to a specific goal and context. Using soil spectroscopy without considering its limitations may lead to unsatisfactory results. Therefore, it is helpful to ask some key questions before designing and implementing a soil spectroscopy project:

| Question | Answer |
|:-----|------:|
| **What soil properties will be predicted?** | Some soil properties are easier to predict than others, depending on the association of spectral features with the soil property of interest. For example, soil organic carbon has a strong correlation with absorbance across the infrared range, while extractable nutrients have a weaker correlation. |
| **What accuracy/precision is required for the project?** | The accuracy and precision of soil spectroscopy estimates depend on the quality and representativeness of the calibration dataset, the spectral range and resolution, and the machine learning method used. Soil spectroscopy may not be suitable for projects that require very high precision, such as detecting small differences in controlled and replicated field trials. However, it may be adequate for projects that aim to classify soils into distinct classes, or to estimate average values from many samples across a defined spatial boundary. |
| **What is the budget for the project?** | Soil spectroscopy can significantly reduce the cost of soil analysis compared to traditional methods, especially when dealing with large numbers of samples. However, costs vary depending on the type and quality of the instrument used. Research-grade bench-top laboratory spectrometers may cost hundreds of thousands of dollars, while lower-cost portable [FTIR](https://en.wikipedia.org/wiki/Fourier-transform_infrared_spectroscopy) instruments may cost only a few thousand. Additional expenses such as sample preparation, data processing, and model development and validation should also be considered. |
| **What instrumentation is available or accessible for the project?** | The choice of instrument depends on the spectral region, resolution, and range most appropriate for the application. Different instruments have different advantages and disadvantages in terms of performance, portability, and usability. For example, VisNIR instruments are more accessible and easier to use than MIR instruments, but may have lower predictive accuracy and precision for some soil properties. |

::: {.callout-tip}
## Tip
By answering these questions, one can assess the feasibility and suitability of soil spectroscopy for a specific project and make informed decisions about the best practices and methods to use.
:::

### Example 1: Cover crop impact on soil organic carbon stocks across Iowa

Cover crops are plants grown in agricultural lands to improve soil health and control erosion. They can also affect soil organic carbon (SOC) stocks, which are important for mitigating climate change and enhancing soil fertility. However, measuring SOC stocks is resource intensive, time-consuming, and expensive, especially at large scales. Soil spectroscopy offers a practical way to assess the impacts of cover crops on SOC stocks — for example, across 20 commercial farms in Iowa.

We would design our study as follows:

- **Requirement**: We need an unbiased estimate of field-level (~20 ha) mean SOC stocks, as well as soil texture components and pH, for each farm.  
- **Study design**: We collect 1 soil core per hectare at 3 depths (0–10 cm, 10–20 cm, and 20–30 cm) for each farm, and divide the farms into two treatments: with cover crops and without cover crops. This results in 2,400 soil samples in total.  
- **Methodology**:
    + We scan all samples with a VisNIR spectrometer to obtain their spectra.  
    + We select 25% of the samples for traditional soil analysis, which involves measuring SOC levels, bulk density, texture, and pH using conventional methods. We use the spectra to subset the most diverse samples, ensuring a representative calibration dataset.  
    + We train a multivariate calibration model using the 25% of samples with reference values and spectra as inputs.  
    + We predict the SOC, texture, and pH of the remaining 75% of samples using only their spectra and the calibration model.

### Example 2: Soil classification using the KSSL MIR Soil Spectral Library

Soil classification is the process of grouping soils into categories based on their physical, chemical, and biological characteristics. It is useful for understanding soil genesis, morphology, distribution, and management. However, obtaining the physicochemical data needed for soil classification is laborious and costly, especially for large numbers of samples. In this example, soil survey has been conducted at state survey offices using compact spectral instruments covering the MIR range. We therefore leverage a large soil spectral library with calibration transfer to obtain the physicochemical characterization data needed for classification.

We may use the following approach:

- We scan all new soil samples with the compact MIR spectrometer to obtain their spectra.  
- The KSSL MIR library already represents samples from soils surveyed across the US.  
- We standardize the instruments using a small set of shared standards to ensure compatibility between the KSSL and the state survey offices' instruments, accounting for small differences in instrument settings and operational procedures. This can be done with routine preprocessing or spectral standardization, depending on the degree of dissimilarity.  
- We assess the representativeness of the KSSL MIR library relative to the new samples using spectral similarity measures. This helps identify whether the new local or regional samples are well-represented by the existing library. If not, localization, transfer learning, or spiking may be necessary.
- We apply a pretrained model from the KSSL MIR to the local sample spectra, with modeling adjustments depending on the previous representation check.  
- We predict several soil properties of interest with uncertainty bounds.  

::: {.callout-tip}
## Tip
If many underrepresented samples are found, we consider sending a portion of them for traditional laboratory analysis and incorporating them into the predictive modeling workflow using localization, transfer learning, or spiking strategies. This helps reduce prediction bias for those new locations.
:::

## Good predictions flow from good data

Predictive soil spectroscopy relies on the quality and consistency of both spectral and reference analytical data. Many factors can affect spectral acquisition or reference analytical quality, and thus impact the performance of predictive models.

From a spectral acquisition standpoint, some of the factors that need to be considered are:

- **Instrument setup and calibration**: Spectral measurements should be made with a well-calibrated and stable instrument, following the manufacturer's instructions and recommendations. Instrument settings, such as spectral range, resolution, and integration time, should be optimized for proper soil characterization and analysis of the soil properties of interest.  
- **Operating environment**: Spectral measurements should be made under controlled and consistent environmental conditions, such as temperature, humidity, lighting, and in some cases under controlled air conditions with chamber purging. Environmental factors can influence the spectral response of soil samples and introduce noise or variability in the data.  
- **Sample preparation**: Soil samples should be prepared in a homogeneous way — including sieving, drying, optional grinding, and surface homogenization. Poor sample preparation can affect the spectral response and bias the estimation of both physical (grain size) and chemical properties (such as organic carbon or mineral occlusion).  
- **Instrument calibration and standard samples**: Spectral instruments should be calibrated prior to use, and measurement sessions should include a consistent reference or standard material — a white reference plate for VisNIR and NIR, or roughened gold or aluminum for MIR — to properly quantify the amount of light reflected by soils. Additional standard soil samples may be included for quality assurance and control to help ensure repeatability and reproducibility of results.

::: {.callout-tip}
## Tip
The soil spectroscopy community has been promoting the use of standard protocols and reference samples for high quality control of laboratory measurements. Examples are provided in the [KSSL Manual Method 7A7](https://www.nrcs.usda.gov/resources/guides-and-instructions/kssl-guidance) and a new international protocol is emerging from the [IEEE SA P4005 working group](https://sagroups.ieee.org/4005/). These reference materials provide guidance for proper spectral acquisition.
:::

## Good practices for model building

Building and evaluating multivariate or machine learning models that relate soil spectra to soil properties of interest requires careful consideration of several factors, such as spectral preprocessing, exploratory data analysis, data splitting, dimensionality reduction, algorithm selection, model interpretation and inspection, and prediction uncertainty estimation. These factors affect the quality, accuracy, and reliability of the models, as well as their interpretability and validity. It is therefore recommended to follow some best practices and assess whether the performance meets the goal and context of the project:

- **Spectral preprocessing**: This step involves removing noise, enhancing the signal, and making the scale comparable across different spectral measurements.

- **Data visualization and exploratory data analysis**: This step aims to answer questions such as: are there any clusters in the data? What factors are driving the spread in the data? Do a few samples have too much influence on the results? Are there any outliers? Data visualization is critical for understanding the structure and patterns of the data.

- **Splitting data for calibration and validation**: This step depends on the goal of the analysis. There are two main approaches: cross-validation and train-test split. For a train-test split, algorithms such as Kennard–Stone or Conditioned Latin Hypercube Sampling can be used to select a representative subset of the data. Cross-validation is suitable when the goal is to use all available data for internal validation and model optimization. A train-test split is more conservative and can provide a better estimate of the model's extrapolation capacity on new data.

- **Dimensionality reduction and feature selection**: This step involves reducing the number of spectral features to speed up computation, avoid redundancy and potential overfitting, and improve interpretability, since spectra are highly dimensional and collinear. Some model types, such as partial least squares regression (PLSR), perform dimensionality reduction during model building, while others do not. Many methods exist for selecting the most important features, such as variable importance in projection (VIP), principal component analysis, or modern variable selection and compression algorithms developed by the machine learning community.

- **Model algorithm choice**: This step is linked to the previous one, as different model types have different assumptions and properties. Many model forms may be applied to spectral data, such as PLSR, Cubist, and neural networks. The choice of model type should be based on experimentation, literature recommendations, performance requirements, and parsimony.

- **Model inspection and interpretation**: This step is essential for evaluating the quality and validity of the model. It involves examining the factors that drive predictions, such as model coefficients and residuals. Some models are easier to interpret than others, depending on their transparency. For example, PLSR can reveal the correlation between spectral features and the soil variable of interest, while deep learning models are more opaque and require additional workarounds for interpretation. Model inspection should also identify any samples that are poorly predicted or show large variation between calibration and validation performance.

- **Outliers and uncertainty**: The final step is to apply the model to new data and assess the accuracy and precision of the predictions. A good practice is to run a small percentage (5–10%) of new samples through the reference analysis method to truly evaluate performance (if not properly done during model splitting), though this may not always be feasible. Alternatively, confidence in the predictions can be increased by: i) checking how well the new samples are represented by the training set using extrapolation metrics such as Hotelling T^2^, inlier distances, or Q statistics; ii) checking whether predictions fall within the range of the measured soil properties; and iii) calculating uncertainty around each prediction point using methods such as [quantile regression](https://en.wikipedia.org/wiki/Quantile_regression) or [conformal prediction](https://en.wikipedia.org/wiki/Conformal_prediction).
