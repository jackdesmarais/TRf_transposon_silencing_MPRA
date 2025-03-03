"""Library analysis module for MPRA data.

This module provides functions and classes for analyzing MPRA library data, including:
- Statistical analysis functions
- Plotting utilities 
- Library class for managing MPRA datasets
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.colors import CenteredNorm

import statsmodels.api as sm
import seaborn as sns
from scipy.stats import pearsonr,gaussian_kde, barnard_exact

from math import ceil
from copy import copy
from bin.utils import PlottingContext

def make_density(x):
    """Calculate kernel density estimate for array.
    
    Args:
        x (array-like): Input array
        
    Returns:
        pandas.Series: Density estimates for input values
    """
    nan_mask = np.logical_not(np.isnan(x))
    x = x[nan_mask]
    if len(x)<2:
        density = pd.Series(np.ones((len(x))))
    else:
        density = pd.Series(gaussian_kde(x).evaluate(x))
    density.index=x.index
    return(density)

def line_plotter(x, y, **kwargs):
    """Plot line through points.
    
    Args:
        x (array-like): x coordinates
        y (array-like): y coordinates
        **kwargs: Additional arguments passed to plt.axline()
    """
    nan_mask = np.logical_not(x.isna()|y.isna())
    if nan_mask.sum()>1:
        x = x[nan_mask]
        y=y[nan_mask]
        ax = plt.gca()
        ax.axline(**kwargs)

def correspondance_plotter(x,y,**kwargs):
    """Plot correlation between x and y values.
    
    Args:
        x (array-like): x coordinates
        y (array-like): y coordinates
        **kwargs: Additional arguments passed to plotting functions
    """
    nan_mask = np.logical_not(x.isna()|y.isna())
    if nan_mask.sum()>1:
        x = x[nan_mask]
        y=y[nan_mask]
        ax = plt.gca()

        r = pearsonr(x,y)
        plt.text(0.1,0.9,f'R:{r.statistic:.2f}',transform=ax.transAxes)

def density_scatter(x,y,**kwargs):
    """Create density scatter plot.
    
    Args:
        x (array-like): x coordinates
        y (array-like): y coordinates
        **kwargs: Additional arguments passed to sns.scatterplot()
    """
    nan_mask = np.logical_not(x.isna()|y.isna())
    if nan_mask.sum()>1:
        x = x[nan_mask]
        y=y[nan_mask]

        if nan_mask.sum()>2:
            d2 = np.stack([x,y])
            density = gaussian_kde(d2).evaluate(d2)
        else:
            density = np.ones(x.shape)

        sort_mask = np.argsort(density)

        sns.scatterplot(x=x.values[sort_mask], y=y.values[sort_mask], hue=density[sort_mask], **kwargs)

def OLS_plotter(x,y,ols_model=None, **kwargs):
    """Plot OLS regression line with confidence intervals.
    
    Args:
        x (array-like): x coordinates
        y (array-like): y coordinates
        ols_model (statsmodels.regression.linear_model.RegressionResultsWrapper, optional): Pre-fit OLS model
        **kwargs: Additional arguments passed to plt.plot()
    """
    nan_mask = np.logical_not(x.isna()|y.isna())
    if nan_mask.sum()>1:
        x = x[nan_mask]
        y=y[nan_mask]

        x_sorter= np.argsort(x.values)
        X = sm.add_constant(x.values[x_sorter])
        
       
        if ols_model is None:
            ols_model = sm.OLS(y.values[x_sorter], X).fit()

        pred = ols_model.get_prediction(X).summary_frame()
        plt.plot(x.values[x_sorter], pred['mean'], **kwargs)
        plt.plot(x.values[x_sorter], pred['mean_ci_lower'], linestyle='--', **kwargs)
        plt.plot(x.values[x_sorter], pred['mean_ci_upper'], linestyle='--', **kwargs)

def gini(sorted_arr):
    """Calculate Gini coefficient for sorted array.
    
    Args:
        sorted_arr (array-like): Sorted input array
        
    Returns:
        float: Gini coefficient
    """
    n = sorted_arr.size
    coef_ = 2. / n
    const_ = (n + 1.) / n
    weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
    return coef_*weighted_sum/(sorted_arr.sum()) - const_

def multi_gini(xs,labels,title='Library skew', scatter_kwargs=dict()):
    """Plot Lorenz curves and calculate Gini coefficients for multiple arrays.
    
    Args:
        xs (list): List of arrays to analyze
        labels (list): Labels for each array
        title (str, optional): Plot title. Defaults to 'Library skew'.
        scatter_kwargs (dict, optional): Additional arguments passed to plt.scatter()
        
    Returns:
        tuple: (list of Gini coefficients, matplotlib figure)
    """
    kwargs = dict(marker='.')
    kwargs.update(scatter_kwargs)
    gs = []
    fig, ax = plt.subplots()
    ax.plot([0,1], [0,1], color='k',label='Perfectly uniform')
    for X_u, label in zip(xs, labels):
        X = X_u.copy()
        X = X[np.logical_not(np.isnan(X))]
        X.sort()
        X_lorenz = X.cumsum() / X.sum()
        X_lorenz = np.insert(X_lorenz, 0, 0)
        gs.append(gini(X))
        ax.scatter(np.arange(X_lorenz.size)/(X_lorenz.size-1), X_lorenz, 
            label=str(label)+' Gini: %.3f'%gs[-1],**kwargs)
    plt.xlabel('Sorted variant fraction')
    plt.ylabel('Cumulative read fraction')
    plt.legend(loc='upper left')
    plt.title(title)
    return(gs, fig)

def depth_plot(xs,labels,title='Library depth'):
    """Plot read depth distributions.
    
    Args:
        xs (list): List of arrays containing read counts
        labels (list): Labels for each array
        title (str, optional): Plot title. Defaults to 'Library depth'.
        
    Returns:
        tuple: (figure with read counts, figure with read fractions)
    """
    fig_read, ax_reads = plt.subplots()
    fig_frac, ax_fraction = plt.subplots()
    for X_u, label in zip(xs, labels):
        sort_obs = sorted(X_u[np.logical_not(np.isnan(X_u))],reverse=True)

        plt.sca(ax_reads)
        plt.loglog(np.arange(len(sort_obs))+1, sort_obs, '.', label=label)

        plt.sca(ax_fraction)
        plt.loglog(np.arange(len(sort_obs))+1, sort_obs/sum(sort_obs), '.', label=label)

    plt.sca(ax_reads)
    plt.xlabel('Variant rank-order')
    plt.ylabel('Variant read count')
    plt.legend(loc='lower left')
    plt.title(title)

    plt.sca(ax_fraction)
    plt.xlabel('Variant rank-order')
    plt.ylabel('Variant fraction')
    plt.legend(loc='lower left')
    plt.title(title)
    return(fig_read, fig_frac)

    

def thresh_collapse(x, thresh=-0.5):
    """Calculate fraction of values above threshold.
    
    Args:
        x (array-like): Input array
        thresh (float, optional): Threshold value. Defaults to -0.5.
        
    Returns:
        float: Fraction of values above threshold
    """
    return((x>thresh).sum()/len(x))

class Library:
    """Class for analyzing MPRA library data.
    
    Attributes:
        data_df (pandas.DataFrame): Library data
        replicates (dict): Replicate information
        id_cols (list): Column names for variant IDs
        group_cols (list): Column names for grouping variants
        rate_method (str): Method for calculating rates
        alphabet (list): Valid characters in sequences
        WT_seq (str): Wild-type sequence
        fitness_measure (str): Column name for fitness values
    """

    def __init__(self, data_df, replicates, id_cols, group_cols, rate_method, alphabet, WT_seq, fitness_measure):
        """Initialize Library object.
        
        Args:
            data_df (pandas.DataFrame): Library data
            replicates (dict): Replicate information
            id_cols (list): Column names for variant IDs
            group_cols (list): Column names for grouping variants
            rate_method (str): Method for calculating rates
            alphabet (list): Valid characters in sequences
            WT_seq (str): Wild-type sequence
            fitness_measure (str): Column name for fitness values
        """
        self.data_df = data_df
        self.replicates = replicates
        self.id_cols = id_cols
        self.group_cols = group_cols
        self.rate_method = rate_method
        self.alphabet = alphabet
        self.WT_seq = WT_seq
        self.fitness_measure = fitness_measure

    @classmethod
    def build_from_reads(cls, data_file, replicates, id_cols, group_cols, rate_method, alphabet, WT_seq, process_call, name_regex, sheet_name=None):
        """Build Library object from read count data file.
        
        Args:
            data_file (str): Path to data file (.xlsx, .csv, or .tsv)
            replicates (dict): Replicate information
            id_cols (list): Column names for variant IDs
            group_cols (list): Column names for grouping variants
            rate_method (str): Method for calculating rates
            alphabet (list): Valid characters in sequences
            WT_seq (str): Wild-type sequence
            process_call (str): Column name for processing
            name_regex (str): Regex pattern for extracting names
            sheet_name (str, optional): Excel sheet name. Defaults to None.
            
        Returns:
            Library: New Library object
            
        Raises:
            ValueError: If data_file is not .xlsx, .csv or .tsv
        """
        if data_file.endswith('.xlsx'):
            data = pd.read_excel(data_file,sheet_name=sheet_name).dropna(axis=0,how='any')
        elif data_file.endswith('.csv'):
            data = pd.read_csv(data_file).dropna(axis=0,how='any')
        elif data_file.endswith('.tsv'):
            data = pd.read_csv(data_file,sep='\t').dropna(axis=0,how='any')
        else:
            raise ValueError('Data file must be an .xlsx, .tsv or .csv file')
        
        normalized_df = pd.DataFrame()
        for id_col in id_cols:
            normalized_df[('meta',id_col)] = data[id_col]

        for replicate, tps in replicates.items():
            for i, tp in enumerate(tps):
                normalized_df[(f'{replicate}_counts',tp)] = data[tp]
                normalized_df[(f'{replicate}_abundance',tp)] = data[tp]/data[tp].sum()
        normalized_df.columns = pd.MultiIndex.from_tuples(normalized_df.columns)

        for replicate, tps in replicates.items():
            if rate_method == 'l2fc':
                normalized_df[('FC',replicate)] = normalized_df[(f'{replicate}_abundance',tps[-1])]/normalized_df[(f'{replicate}_abundance',tps[0])]
                normalized_df[('l2fc',replicate)] = np.log2(normalized_df[('FC',replicate)])
            elif rate_method == 'l10fc':
                normalized_df[('FC',replicate)] = normalized_df[(f'{replicate}_abundance',tps[-1])]/normalized_df[(f'{replicate}_abundance',tps[0])]
                normalized_df[('l10fc',replicate)] = np.log10(normalized_df[('FC',replicate)])
        normalized_df = normalized_df.sort_index(axis=1)
        meta = normalized_df[('meta',process_call)].str.extract(name_regex)
        for col in meta.columns:
            normalized_df[('meta',col)] = meta[col]
        normalized_df = normalized_df.sort_index(axis=1)

        self = cls(normalized_df, replicates, id_cols, group_cols, rate_method, alphabet, WT_seq, rate_method)
        return(self)
    
    def save(self, filename):
        """Save Library object to a pickle file.
        
        Args:
            filename (str): Path to save the library object
            
        Raises:
            ValueError: If filename does not end in .pkl
        """
        if not filename.endswith('.pkl'):
            raise ValueError('Filename must end in .pkl')
            
        import pickle
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
    
    @classmethod
    def load(cls, filename):
        """Load Library object from a pickle file.
        
        Args:
            filename (str): Path to the library pickle file
            
        Returns:
            Library: Loaded Library object
            
        Raises:
            ValueError: If filename does not end in .pkl
            ValueError: If loaded object is not a Library instance
        """
        if not filename.endswith('.pkl'):
            raise ValueError('Filename must end in .pkl')
            
        import pickle
        with open(filename, 'rb') as f:
            loaded = pickle.load(f)
            
        if not isinstance(loaded, cls):
            raise ValueError(f'Loaded object is not a {cls.__name__} instance')
            
        return loaded
    
    
    @property
    def fitness_frame(self):
        """Get fitness data frame.
        
        Returns:
            pandas.DataFrame: Fitness values for each variant
        """
        df = self.data_df.set_index([('meta',col) for col in self.group_cols+self.id_cols])[self.fitness_measure].rename(mapper=lambda x: f'{self.fitness_measure} {x}', 
                axis='columns',
                level=0,
                inplace=False)
        df.index.names = self.group_cols+self.id_cols
        return(df.reset_index())
    
    @property
    def rep_average_frame(self):
        """Get replicate-averaged fitness data frame.
        
        Returns:
            pandas.DataFrame: Mean and std of fitness values across replicates
        """
        df = self.fitness_frame
        df = df.set_index(self.group_cols+self.id_cols).agg(['mean','std'],axis=1).rename(mapper=lambda x: f'{self.fitness_measure} {x}', 
                axis='columns',
                level=0,
                inplace=False)
        return(df.reset_index())
    
    @property
    def group_average_frame(self):
        """Get group-averaged fitness data frame.
        
        Returns:
            pandas.DataFrame: Mean and std of fitness values across groups
        """
        df = self.fitness_frame
        df = df.groupby(self.group_cols).agg(['mean','std'],axis=1).rename(mapper=lambda x: f'{self.fitness_measure} {x}', 
                axis='columns',
                level=1,
                inplace=False)
        return(df.reset_index())

    @property
    def total_average_frame(self):
        """Get total averaged fitness data frame.
        
        Returns:
            pandas.DataFrame: Mean fitness values across all replicates and groups
        """
        df = self.rep_average_frame
        df = df.groupby(self.group_cols)[f'{self.fitness_measure} mean'].mean().rename(f'{self.fitness_measure} mean')
        return(df.reset_index())
    
    def check_skew(self, sns_settings='notebook', rc_settings={}, title='Library skew', palette=None):
        """Check library skew using Lorenz curves.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            title (str, optional): Plot title. Defaults to 'Library skew'.
            
        Returns:
            matplotlib.figure.Figure: Lorenz curve plot
        """
        with PlottingContext(sns_settings=sns_settings, rc_settings=rc_settings, palette=palette):

            ginis, fig = multi_gini([self.data_df[(f'{key}_counts', tp)][~self.data_df[(f'{key}_counts', tp)].isna()].values for key, tps in self.replicates.items() for tp in tps],
                             [f'{key} {tp}' for key, tps in self.replicates.items() for tp in tps],
                             title=title)
            
            
        return(fig)
    
    def check_depth(self, sns_settings='notebook', rc_settings={}, title='Read depth', palette=None):
        """Check read depth distributions.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            title (str, optional): Plot title. Defaults to 'Read depth'.
            
        Returns:
            tuple: (figure with read counts, figure with read fractions)
        """
        with PlottingContext(sns_settings=sns_settings, rc_settings=rc_settings, palette=palette):

            fig_read, fig_frac = depth_plot([self.data_df[(f'{key}_counts', tp)][~self.data_df[(f'{key}_counts', tp)].isna()].values for key, tps in self.replicates.items() for tp in tps],
                             [f'{key} {tp}' for key, tps in self.replicates.items() for tp in tps],
                             title=title)
            
        return(fig_read, fig_frac)
    
    def check_barcode_distribution_skew(self,sns_settings='notebook', threshold=0, rc_settings={}, title='Read depth', palette=None):
        """Check skew in barcode distribution using Lorenz curves.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            threshold (int, optional): Minimum read count threshold. Defaults to 0.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            title (str, optional): Plot title. Defaults to 'Read depth'.
            
        Returns:
            matplotlib.figure.Figure: Lorenz curve plot
        """
        with PlottingContext(sns_settings=sns_settings, rc_settings=rc_settings, palette=palette):
            barcode_counts_df = self.data_df.groupby([('meta',col) for col in self.group_cols])[[(k+'_counts', col) for k, v in self.replicates.items() for col in v]].apply(count_over_thresh, threshold=threshold)
            ginis, fig = multi_gini([barcode_counts_df[(k+'_counts', col)][~barcode_counts_df[(f'{k}_counts', col)].isna()].values for k, v in self.replicates.items() for col in v],
                             [f'{key} {tp}' for key, tps in self.replicates.items() for tp in tps],
                             title=title)
            
        return(fig)

    def check_barcode_distribution_depth(self, sns_settings='notebook', threshold=0, rc_settings={}, title='Read depth', palette=None):
        """Check read depth distributions for barcodes.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            threshold (int, optional): Minimum read count threshold. Defaults to 0.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            title (str, optional): Plot title. Defaults to 'Read depth'.
            
        Returns:
            tuple: (figure with read counts, figure with read fractions)
        """
        with PlottingContext(sns_settings=sns_settings, rc_settings=rc_settings, palette=palette):
            barcode_counts_df = self.data_df.groupby([('meta',col) for col in self.group_cols])[[(k+'_counts', col) for k, v in self.replicates.items() for col in v]].apply(count_over_thresh, threshold=threshold)
            fig_read, fig_frac = depth_plot([barcode_counts_df[(k+'_counts', col)][~barcode_counts_df[(f'{k}_counts', col)].isna()].values for k, v in self.replicates.items() for col in v],
                             [f'{key} {tp}' for key, tps in self.replicates.items() for tp in tps],
                             title=title)
            
        return(fig_read, fig_frac)
    
    
    def _make_replicate_checks(self, to_plot, sns_settings='notebook', rc_settings={}, grid_kwargs=dict(), palette=None):
        """Generate replicate-level fitness plots.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            grid_kwargs (dict, optional): Additional arguments for the grid. Defaults to {}.    
            palette (str, optional): Palette to use. Defaults to None.
            
        Returns:
            tuple: (unaveraged plot, averaged plot)
        """
        with PlottingContext(sns_settings=sns_settings, rc_settings=rc_settings, palette=palette):
            if len(to_plot.columns)==2:
                to_plot = to_plot.dropna(how='any',axis='index')
                xname = list(to_plot.columns)[0]
                yname = list(to_plot.columns)[1]

                unaveraged_g = sns.JointGrid(data = to_plot, x=xname, y=yname, **grid_kwargs) 
                unaveraged_g.plot_joint(density_scatter,edgecolor="none", legend=False)
                unaveraged_g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
                unaveraged_g.plot_joint(correspondance_plotter)
                unaveraged_g.plot_marginals(sns.histplot)   
                unaveraged_g.plot_joint(OLS_plotter, color='C0', alpha=0.3)
            elif len(to_plot.columns)>2:
                unaveraged_g = sns.PairGrid(to_plot, corner=True, **grid_kwargs)
                unaveraged_g.map_diag(sns.histplot)
                unaveraged_g.map_lower(density_scatter,edgecolor="none", legend=False)
                unaveraged_g.map_lower(correspondance_plotter)
                unaveraged_g.map_lower(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
                unaveraged_g.map_lower(OLS_plotter, color='C0', alpha=0.3)
            else:
                raise ValueError('to_plot must have at least 2 columns')
                
            
        return(unaveraged_g)
            
    
    def make_replicate_checks(self, sns_settings='notebook', rc_settings={}, grid_kwargs=dict(), palette=None):
        """Generate replicate-level fitness plots.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            grid_kwargs (dict, optional): Additional arguments for the grid. Defaults to {}.
            
        Returns:
            tuple: (unaveraged plot, averaged plot)
        """
        return(self._make_replicate_checks(self.data_df[self.fitness_measure], sns_settings, rc_settings, grid_kwargs, palette))
                

    def make_grouped_replicate_plots(self, grouping_dict=None, sns_settings='notebook', rc_settings={}, grid_kwargs=dict(), palette=None):
        """Generate grouped replicate-level fitness plots.
        
        Args:
            grid_kwargs (dict): Additional arguments for the grid
            
        Returns:
            seaborn plot: Either JointGrid or PairGrid depending on number of replicates
        """

        assert grouping_dict is not None, 'grouping_dict must be provided'

        rep_average_ff = pd.concat([self.data_df[self.fitness_measure][to_pool].mean(axis=1).rename(name).to_frame() for name, to_pool in grouping_dict.items()],axis='columns')
        return(self._make_replicate_checks(rep_average_ff, sns_settings, rc_settings, grid_kwargs, palette))



    def make_initial_skew_checks(self, sns_settings='notebook', rc_settings={}, palette=None):
        """Make plots checking initial abundance vs fitness measure skew.
        
        Args:
            sns_settings (str, optional): Seaborn plotting settings. Defaults to 'notebook'.
            rc_settings (dict, optional): Dictionary of rc parameters to use. Defaults to {}.
            
        Returns:
            tuple: A tuple containing:
                - self (Library): The Library object
                - figs (list): List of seaborn.JointGrid plots showing skew
        """
        with PlottingContext(sns_settings=sns_settings, rc_settings=rc_settings, palette=palette):
            figs = []
            for key, tps in self.replicates.items():
                to_plot = self.data_df[[(f'{key}_abundance', tps[0]),(self.fitness_measure, key)]].dropna(how='any',axis='index')
                to_plot.columns = [f'{tps[0]} abundance', f'{key} {self.fitness_measure}']

                g = sns.JointGrid(data = to_plot, x=f'{tps[0]} abundance', y=f'{key} {self.fitness_measure}') 
                g.plot_joint(density_scatter, edgecolor="none", legend=False)
                # g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
                g.plot_joint(correspondance_plotter)
                g.plot_marginals(sns.histplot) 
                figs.append(g)

        return(self, figs)