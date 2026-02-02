
## Overview

In optimisation cases, we face many problems: correlation, structural breaks, outliers; the hypotheses of a good linear model are rarely verified in real-life data. That is where the role of different algorithms appears.

Statistical methods like Stepwise forward and backward, and machine learning algorithms like LASSO, LARS, and Elastic Net are used in variable selection, but do we really know the performance of these models, and in which cases we should use one algorithm rather than another?

This study aims to answer these questions by evaluating these algorithms under different data scenarios.

We test all these algorithms in several cases using Monte Carlo simulations implemented in SAS. For each algorithm and for each scenario (outliers, internal correlation, independence, structural breaks), we evaluate model selection performance and calculate probabilities of overfitting, underfitting, and perfect model selection.

## Data generation

First, we defined several functions to generate synthetic data under different scenarios and store them in dedicated tables within a library. These scenarios include independence, structural breaks, outliers, and collinearity.

The data generation process follows the methodology used in the course, where different data-generating processes (DGPs) simulate realistic violations of classical linear model assumptions.

This allows us to test algorithm robustness under controlled but realistic conditions.

## Variable selection methods

In the second step, we define functions implementing different selection methods:

* Stepwise selection (forward and backward),
* Stagewise regression,
* LASSO (with lambda chosen by cross-validation or fixed manually),
* LARS,
* Elastic Net.

Each method is applied repeatedly on the generated datasets through Monte Carlo simulations.

## Performance evaluation

At each simulation iteration, several indicators are computed:

* Overfitting: at least one inactive variable is selected.
* Underfitting: at least one active variable is not selected.
* Perfect fitting: all relevant variables are correctly selected.
* Training MSE and Test MSE.
* Overfitting ratio measured by MSE_train / MSE_test.
* False positive and false negative rates.
* Bias and variance of estimated coefficients.

## Analysis and comparison

Results are summarized through comparative tables and visualizations, including bar plots comparing bias and variance across methods and scenarios.

The objective is to compare the performance and stability of each method under different data conditions and provide guidance on which algorithm performs best depending on the scenario.

## Goal of the project

The final goal is to better understand when each variable selection algorithm should be used and how robust these methods are when classical modeling assumptions are violated.

