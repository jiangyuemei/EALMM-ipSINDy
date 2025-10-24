# EAAM-ipSINDy: Ensemble Adaptive Adams-Moulton with Inner-Product Sparse Identification

A robust framework for discovering nonlinear dynamical systems from noisy data using ensemble adaptive libraries and Adams-Moulton linear multistep methods.

## 🌟 Overview

EAAM-ipSINDy implements a novel approach for sparse identification of nonlinear dynamics that combines:

- **Adams-Moulton linear multistep methods** for robust derivative approximation
- **Ensemble adaptive libraries** for parsimonious feature selection  
- **Inner-product driven sparse regression** for noise-resilient coefficient identification
- **Adaptive moving average filtering** for effective noise reduction

This method significantly improves upon traditional SINDy and WSINDy approaches in both accuracy and computational efficiency, particularly under high-noise conditions.

## 🚀 Features

- **Robust to Noise**: Handles noise levels up to 50% signal-to-noise ratio
- **Multi-order Support**: Adams-Moulton methods of orders 1-5
- **Ensemble Processing**: 200+ ensemble realizations for statistical reliability
- **Adaptive Libraries**: Dynamically constructs optimal feature sets
- **Cross-validation**: Ensures generalization and prevents overfitting

## 📋 Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox
- Optimization Toolbox (recommended)

## 🛠 Installation

1. Clone the repository:
```bash
git clone https://github.com/jiangyuemei/EAAM-ipSINDy.git
cd EAAM-ipSINDy
