# 🧬 Biological Random Walks (BRW)

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Streamlit](https://img.shields.io/badge/Streamlit-Web%20App-red.svg)](https://streamlit.io/)

> **Integrating Gene Expression for Tissue-Specific Prediction**

A comprehensive framework for disease gene prioritization using Biological Random Walks algorithm with interactive GUI and batch processing capabilities.

## 🚀 Quick Start

### Prerequisites
- Python 3.7 or higher
- pip package manager

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/BiologicalRandomWalks.git
   cd BiologicalRandomWalks
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Download required data files**
   - Place your data files in the `data_set/` directory
   - Ensure `Dataset OncoKB.xlsx` is in the root directory

## 🖥️ Usage Options

### Option 1: Interactive GUI (Recommended for beginners)
```bash
streamlit run app_gui.py
```
- User-friendly web interface
- Interactive parameter configuration
- Real-time result visualization
- No command-line knowledge required

### Option 2: Command Line
```bash
python main.py -p <ppi_network> -s <seed_set> -a <ontology> -do <disease_ontology> -o <output>
```

### Option 3: Batch Processing
```bash
python run_brw_batch.py
```
- Automated processing of all cancer types
- Comprehensive parameter testing
- Generates summary reports

## 📊 Features

- **Multiple Data Integration**: PPI networks, co-expression, DE genes, ontologies
- **Flexible Configuration**: Ablation studies and parameter tuning
- **Cancer-Specific Analysis**: Support for 7 cancer types (BRCA, COAD, LUAD, THCA, BLCA, PRAD, STAD)
- **OncoKB Validation**: Automatic gene prioritization validation
- **Result Visualization**: Charts, heatmaps, and statistical analysis
- **Export Functionality**: CSV downloads and publication-ready figures

## 🏗️ Project Structure

```
BiologicalRandomWalks/
├── app_gui.py                 # Streamlit GUI application
├── main.py                    # Main BRW execution script
├── run_brw_batch.py          # Batch processing script
├── check_genes.py            # OncoKB validation
├── biological_random_walks/   # Core BRW algorithm
├── data_preprocessing/        # Data preparation scripts
├── data_set/                  # Input data directory
├── outputs/                   # Results directory
├── requirements.txt           # Python dependencies
└── README.md                  # Detailed documentation
```

## 🔬 Algorithm Overview

Biological Random Walks (BRW) extends the Random Walk with Restart algorithm by:

1. **Computing disease-specific annotations** from known disease genes
2. **Individual teleporting probabilities** based on biological similarity
3. **Weighted PPI interactions** using node annotations
4. **Biased flow propagation** toward functionally related genes

## 📚 Documentation

- **Full Documentation**: See [README.md](README.md) for comprehensive details
- **GUI Guide**: Interactive web interface documentation
- **API Reference**: Command-line options and parameters
- **Examples**: Sample data and usage scenarios

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Original Authors**: Michele Gentili, Leonardo Martini, Manuela Petti, Luca Becchetti, Lorenzo Farina, Marialuisa Sponziello
- **Institution**: Department of Computer, Control, and Management Engineering Antonio Ruberti, Sapienza University of Rome
- **Research**: Translational and Precision Medicine Department, Sapienza University of Rome

## 📞 Contact

- **Project Issues**: [GitHub Issues](https://github.com/yourusername/BiologicalRandomWalks/issues)
- **Original Authors**: 
  - Michele Gentili: gentili@diag.uniroma1.it
  - Leonardo Martini: martini@diag.uniroma1.it

---

**⭐ Star this repository if you find it useful!**

**Happy Analyzing with BRW! 🧬✨**
