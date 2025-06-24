import plotly.graph_objects as go
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
from helper import *

def create_html_file(plot_html_list, output_file="combined_plots.html"):
    """
    Combines multiple Plotly plots into a single HTML file.

    Parameters:
    - plot_html_list: List of HTML content for each plot.
    - output_file: Name of the output HTML file.
    """
    # Start the HTML content
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Combined Plots</title>
        <style>
            .plot-container {
                margin-bottom: 40px;
            }
        </style>
    </head>
    <body>
        <h1>Combined Plots</h1>
    """

    # Add each plot to the HTML content
    for plot_html in plot_html_list:
        html_content += f"""
        <div class="plot-container">
            {plot_html}
        </div>
        """

    # Close the HTML content
    html_content += """
    </body>
    </html>
    """

    # Write the HTML content to a file
    with open(output_file, "w") as f:
        f.write(html_content)

    print(f"HTML file saved as {output_file}")


# Function to create an interactive plot and save it as HTML
def plot_CLL_count(patient_id,times_sample, CLL_count, UMI_start, UMI_end, treatment_start, treatment_end, ):
    # Create the plot
    fig = go.Figure()

    # Add the main CLL count line
    fig.add_trace(go.Scatter(
        x=times_sample,
        y=CLL_count,
        mode='lines+markers',
        name='CLL Count',
        line=dict(color='blue')
    ))

    # Highlight UMI timepoints
    fig.add_trace(go.Scatter(
        x=times_sample[UMI_start:UMI_end],
        y=CLL_count[UMI_start:UMI_end],
        mode='markers',
        name='UMI Timepoints',
        marker=dict(color='red', size=10, symbol='circle')
    ))

    # Add a shaded region for the treatment period
    fig.add_vrect(
        x0=treatment_start,
        x1=treatment_end,
        fillcolor="lightgray",
        opacity=0.5,
        layer="below",
        line_width=0,
        annotation_text="Fludarabine Treatment",
        annotation_position="top left"
    )

    # Customize the layout
    fig.update_layout(
        title=f'CLL Count Over Time for Patient {patient_id}',
        xaxis_title='Days from Diagnosis',
        yaxis_title='CLL Count (log scale)',
        yaxis_type="log",
        showlegend=True,
        template="plotly_white"
    )

    # Return html content
    return fig.to_html(full_html=False)

def plot_metadata_table(df, patient_id):
    # Create a Plotly table
    table = go.Table(
        header=dict(
            values=list(df.columns),  # Column names
            fill_color='lightblue',   # Header background color
            align='left'              # Align text to the left
        ),
        cells=dict(
            values=[df[col] for col in df.columns],  # Cell values
            fill_color='white',                      # Cell background color
            align='left'                             # Align text to the left
        )
    )

    # Create a figure and add the table
    fig = go.Figure(data=[table])

    # Update layout for better appearance
    fig.update_layout(
        title=f'Patient ID:{patient_id}',
        title_x=0.5,  # Center the title
        margin=dict(l=20, r=20, t=40, b=20)  # Adjust margins
    )

    return fig.to_html(full_html=False)


def plot_tree_plotly(tree_df, tree_selected, edge_labels=None):
    edges = tree_df.loc[tree_selected, 'edges'].split(',')
    cluster_list = []
    for i in edges:
        new_list = i.split('-')
        for j in new_list:
            if (j != 'None') & (j not in cluster_list):
                cluster_list.append(j)
    cluster_list = [int(i) for i in cluster_list]

    DG = nx.DiGraph()
    for edge in edges:
        nodes = edge.split('-')
        if nodes[0] != 'None':
            DG.add_edge(int(nodes[0]), int(nodes[1]))

    pos = graphviz_layout(DG, prog='dot')

    # Create edge traces
    edge_traces = []
    for edge in DG.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_color = ClusterColors.get_hex_string(edge[1])

        edge_trace = go.Scatter(
            x=[x0, x1, None],
            y=[y0, y1, None],
            line=dict(width=2, color=edge_color),
            hoverinfo='none',
            mode='lines')
        edge_traces.append(edge_trace)

    # Create node trace
    node_x = []
    node_y = []
    node_text = []
    node_color = []
    for node in DG.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(str(node))
        node_color.append(ClusterColors.get_hex_string(node))

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=node_text,
        textposition="middle center",
        textfont=dict(color='white', size=14),
        hoverinfo='text',
        marker=dict(
            color=node_color,
            size=40,
            line=dict(width=2, color='DarkSlateGrey'))
    )

    # Create figure
    fig = go.Figure(data=edge_traces + [node_trace],
                    layout=go.Layout(
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    return fig.to_html(full_html=False)


def plot_ccf_tree_combined(tree_df, tree_selected, ccf_df, times_sample, treatment_df):
    fig = make_subplots(
        rows=1,
        cols=2,
        column_widths=[0.5, 0.5],
        horizontal_spacing=0.05,  # Minimal gap
        subplot_titles=("Phylogenetic Tree", "Cancer Cell Fraction (CCF)")
    )

    # --- Tree Plot (Left Subplot) ---
    edges = tree_df.loc[tree_selected, 'edges'].split(',')
    cluster_list = []
    for i in edges:
        new_list = i.split('-')
        for j in new_list:
            if (j != 'None') & (j not in cluster_list):
                cluster_list.append(j)
    cluster_list = [int(i) for i in cluster_list]

    DG = nx.DiGraph()
    for edge in edges:
        nodes = edge.split('-')
        if nodes[0] != 'None':
            DG.add_edge(int(nodes[0]), int(nodes[1]))

    pos = graphviz_layout(DG, prog='dot')

    # Edge traces
    for edge in DG.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_color = ClusterColors.get_hex_string(edge[1])
        fig.add_trace(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                line=dict(width=2, color=edge_color),
                hoverinfo='none',
                mode='lines',
                showlegend=False
            ),
            row=1, col=1
        )

    # Node trace
    node_x = []
    node_y = []
    node_text = []
    node_color = []
    for node in DG.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(str(node))
        node_color.append(ClusterColors.get_hex_string(node))

    fig.add_trace(
        go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            text=node_text,
            textposition="middle center",
            textfont=dict(color='white', size=14),
            hoverinfo='text',
            marker=dict(
                color=node_color,
                size=40,
                line=dict(width=2, color='DarkSlateGrey')
            ),
            showlegend=False
        ),
        row=1, col=1
    )

    # --- 2. CCF Plot (Right Subplot) ---
    cols = ['Sample_ID', 'Cluster_ID', 'postDP_ccf_mean', 'postDP_ccf_CI_low', 'postDP_ccf_CI_high']
    ccf_df = ccf_df[cols]
    cluster_list_ccf = ccf_df.Cluster_ID.unique().tolist()
    x_axis = [i / 365 for i in times_sample]  # Days to years

    for cluster in cluster_list_ccf:
        cluster_data = ccf_df[ccf_df.Cluster_ID == cluster]
        y_mean = cluster_data.postDP_ccf_mean
        y_ci_low = cluster_data.postDP_ccf_CI_low
        y_ci_high = cluster_data.postDP_ccf_CI_high

        # Create the confidence interval as a single trace
        # Combine x values (forward and backward) and y values (high and low)
        x_combined = list(x_axis) + list(x_axis[::-1])  # x values forward then backward
        y_combined = list(y_ci_high) + list(y_ci_low[::-1])  # high values then low values reversed

        # Plot confidence interval as filled area
        fig.add_trace(
            go.Scatter(
                x=x_combined,
                y=y_combined,
                fill='toself',
                fillcolor=f'rgba({",".join(str(int(ClusterColors.get_hex_string(cluster)[i:i + 2], 16)) for i in (1, 3, 5))}, 0.2)',
                # Convert hex to rgba with transparency
                line=dict(color='rgba(255,255,255,0)'),  # Invisible line
                showlegend=False,
                legendgroup=f'Cluster {cluster}',
                hoverinfo='skip'
            ),
            row=1, col=2
        )

        # Mean CCF
        fig.add_trace(
            go.Scatter(
                x=x_axis,
                y=y_mean,
                mode='lines+markers',
                line=dict(color=ClusterColors.get_hex_string(cluster)),
                marker=dict(color=ClusterColors.get_hex_string(cluster)),
                name=f'Cluster {cluster}',
                legendgroup=f'Cluster {cluster}'
            ),
            row=1, col=2
        )

    # Treatment annotations
    for _, row_data in treatment_df.iterrows():
        treatment_name = row_data.tx
        start = row_data.tx_start / 365
        end = row_data.tx_end / 365 if not np.isnan(row_data.tx_end) else max(x_axis)
        fig.add_vrect(
            x0=start, x1=end,
            fillcolor="lightgray", opacity=0.2,
            layer="below", line_width=0,
            annotation_text=treatment_name,
            annotation_position="top left",
            row=1, col=2
        )

    # --- Layout Adjustments ---
    # Tree subplot: Hide axes (matches nx.draw)
    fig.update_xaxes(showgrid=False, showticklabels=False, row=1, col=1)
    fig.update_yaxes(showgrid=False, showticklabels=False, row=1, col=1)

    # CCF subplot: Add grid/labels
    fig.update_xaxes(title_text="Time (years)", showgrid=True, gridcolor='lightgray', row=1, col=2)
    fig.update_yaxes(title_text="CCF", showgrid=True, gridcolor='lightgray', row=1, col=2)

    fig.update_layout(
        height=600,
        margin=dict(l=20, r=20, b=20, t=40),
        legend=dict(x=1.02, y=1)  # Legend outside
    )

    # Save to HTML
    return fig.to_html(full_html=False)


def plot_subclones(clusters, times_sample, CLL_count_sample, log_subclone_sample, extrapolate_start_idx, times_aft_tx,
                   treatment_df, treatment_end, ):
    """
    Plot subclones and extrapolate their behavior after treatment using Plotly.

    Args:
        clusters (list): List of cluster IDs to plot.
        times_sample (list): List of timepoints since treatment start.
        CLL_count (list): List of CLL count estimates.
        log_subclone_sample (dict): Dictionary of log subclone counts for each cluster.
        extrapolate_start_idx (int): Index to start extrapolation.
        times_aft_tx (list): Timepoints after treatment for extrapolation.
        treatment (pd.DataFrame): Treatment data.
        ClusterColors: Object to get cluster colors.
    """
    # Create a subplot with two rows
    fig = go.Figure()

    # Plot total WBC
    x_year = [i / 365 for i in np.array(times_sample)]
    fig.add_trace(
        go.Scatter(x=x_year, y=np.log(CLL_count_sample), mode='markers', name='CLL count estimate', marker=dict(color='red')),

    )

    # Plot subclones and extrapolate after treatment
    predicted_end_tx_clones = []
    for i in clusters:
        y_sub = np.array(log_subclone_sample[i])
        fig.add_trace(
            go.Scatter(x=x_year, y=y_sub, mode='markers', name=f'Cluster {i}',
                       marker=dict(color=ClusterColors.get_hex_string(i))),

        )

        # Extrapolate after treatment
        extrapolate_times = times_sample[extrapolate_start_idx:]
        extrapolate_year = [i / 365 for i in np.array(extrapolate_times)]
        extrapolate_subclone = log_subclone_sample[i][extrapolate_start_idx:]

        # Fit a linear model for extrapolation
        linear_model = np.polyfit(extrapolate_year, extrapolate_subclone, 1)
        predict_year = [i / 365 for i in np.array(times_aft_tx)]
        predicted = np.polyval(linear_model, predict_year)
        predicted_end_tx_clones.append(predicted[0])

        # Plot the extrapolated line
        fig.add_trace(
            go.Scatter(x=predict_year, y=predicted, mode='lines', name=f'Cluster {i} (extrapolated)',
                       line=dict(color=ClusterColors.get_hex_string(i))),

        )

    # Plot treatment effects
    times_during_tx = [0, treatment_end / 365]
    tx_start_clones = [log_subclone_sample[i][0] for i in clusters]

    for i in clusters:
        fig.add_trace(
            go.Scatter(x=times_during_tx, y=[tx_start_clones[i - 1], predicted_end_tx_clones[i - 1]],
                       mode='lines+markers', name=f'Cluster {i} (treatment)',
                       line=dict(color=ClusterColors.get_hex_string(i))),

        )

    # Add treatment annotations
    for i, row in treatment_df.iterrows():
        treatment_name = row.tx
        start = row.tx_start / 365
        end = treatment_end / 365
        fig.add_vrect(x0=start, x1=end, fillcolor="lightgray", opacity=0.5, line_width=0,
                      annotation_text=treatment_name, annotation_position="top left", )

    # Update layout
    fig.update_layout(
        title="CLL Count and Subclones with Treatment Effects",
        xaxis_title="Time (years)",
        yaxis_title="log CLL count estimates",

        showlegend=True,
        height=800,
    )

    # Return HTML content
    return fig.to_html(full_html=False)


def plot_linear_model_mcmc(clusters, times_sample, CLL_count_sample, log_subclone_sample_mcmc_with_uniform_noise, extrapolate_start_idx, times_aft_tx, treatment_df, treatment_end):
    """
    Plots subclones with MCMC data for all clusters using Plotly.
    Includes both growth rate and decay rate histograms and tables.

    Parameters:
    - clusters: List of cluster IDs to plot.
    """
    num_clusters = len(clusters)
    # Calculate rows needed:
    # - One row per pair of clusters
    # - One row for histograms
    # - Two rows for tables (one for growth, one for decay)
    num_rows = ((num_clusters + 1) // 2) + 3  # +3 for histograms and tables
    num_cols = 2

    # Define subplot specs
    specs = []
    # Add rows for cluster plots
    for _ in range((num_clusters + 1) // 2):
        specs.append([{'type': 'xy'}, {'type': 'xy'}])
    # Add row for histograms
    specs.append([{'type': 'xy'}, {'type': 'xy'}])
    # Add rows for tables
    specs.append([{'type': 'table', 'colspan': 2}, None])  # Growth table
    specs.append([{'type': 'table', 'colspan': 2}, None])  # Decay table

    # Create subplots with appropriate titles
    subplot_titles = [f'Cluster {cluster}' for cluster in clusters]
    subplot_titles.extend([""] if len(clusters) % 2 == 1 else [])
    subplot_titles.extend(['Growth Rate Histogram', 'Decay Rate Histogram'])

    fig = make_subplots(
        rows=num_rows, cols=num_cols,
        subplot_titles=subplot_titles,
        vertical_spacing=0.1,
        specs=specs,
        row_heights=[1] * ((num_clusters + 1) // 2) + [1, 0.5, 0.5]  # Adjust table row heights
    )

    # Plot subclones for each cluster
    slopes_mcmc = {}
    slopes_mcmc_decay = {}
    for idx, cluster in enumerate(clusters):
        row = (idx // 2) + 1
        col = (idx % 2) + 1

        # Plot total WBC
        x_year = [i / 365 for i in times_sample]
        fig.add_trace(
            go.Scatter(
                x=x_year,
                y=np.log(CLL_count_sample),
                mode='markers',
                marker=dict(color='red'),
                name='CLL count estimate'
            ),
            row=row, col=col
        )

        # Plot subclones and extrapolate
        slopes_mcmc[cluster] = []
        slopes_mcmc_decay[cluster] = []
        for iter_idx in range(250):
            y_sub = log_subclone_sample_mcmc_with_uniform_noise[cluster][iter_idx]
            fig.add_trace(
                go.Scatter(
                    x=x_year,
                    y=y_sub,
                    mode='markers',
                    marker=dict(color=ClusterColors.get_hex_string(cluster)),
                    showlegend=False
                ),
                row=row, col=col
            )

            # Extrapolate after treatment
            extrapolate_times = times_sample[extrapolate_start_idx:]
            extrapolate_year = [i / 365 for i in extrapolate_times]
            extrapolate_subclone = y_sub[extrapolate_start_idx:]
            linear_model = np.polyfit(extrapolate_year, extrapolate_subclone, 1)
            slopes_mcmc[cluster].append(linear_model[0])

            predict_year = [i / 365 for i in times_aft_tx]
            predicted = np.polyval(linear_model, predict_year)

            # Calculate decay rate
            slope_decay = (predicted[0] - y_sub[0]) / (treatment_end / 365)
            slopes_mcmc_decay[cluster].append(slope_decay)

            fig.add_trace(
                go.Scatter(
                    x=predict_year,
                    y=predicted,
                    mode='lines',
                    line=dict(color=ClusterColors.get_hex_string(cluster)),
                    showlegend=False
                ),
                row=row, col=col
            )

            # Plot during treatment
            times_during_tx = [0, treatment_end / 365]
            tx_start_clones = y_sub[0]
            fig.add_trace(
                go.Scatter(
                    x=times_during_tx,
                    y=[tx_start_clones, predicted[0]],
                    mode='lines+markers',
                    line=dict(color=ClusterColors.get_hex_string(cluster)),
                    marker=dict(color=ClusterColors.get_hex_string(cluster)),
                    showlegend=False
                ),
                row=row, col=col
            )

        # Add treatment information
        for i, row_data in treatment_df.iterrows():
            treatment_name = row_data.tx
            start = row_data.tx_start / 365
            end = treatment_end / 365 if not np.isnan(row_data.tx_end) else max(x_year)
            fig.add_vrect(
                x0=start, x1=end,
                fillcolor="gray", opacity=0.5,
                layer="below", line_width=0,
                annotation_text=treatment_name, annotation_position="top left",
                row=row, col=col
            )

        # Set axis labels
        fig.update_xaxes(title_text="Time (years)", row=row, col=col)
        fig.update_yaxes(title_text="log CLL count estimates (×10^6 cells/ml)", row=row, col=col)
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', row=row, col=col)
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', row=row, col=col)

    # Plot growth rate histogram
    row_hist = num_rows - 2
    for cluster in clusters:
        fig.add_trace(
            go.Histogram(
                x=slopes_mcmc[cluster],
                name=f'Cluster {cluster} Growth',
                marker_color=ClusterColors.get_hex_string(cluster),
                opacity=0.5
            ),
            row=row_hist, col=1
        )
    fig.update_xaxes(title_text="Growth Rate (slope)", row=row_hist, col=1)
    fig.update_yaxes(title_text="Count", row=row_hist, col=1)

    # Plot decay rate histogram
    for cluster in clusters:
        fig.add_trace(
            go.Histogram(
                x=slopes_mcmc_decay[cluster],
                name=f'Cluster {cluster} Decay',
                marker_color=ClusterColors.get_hex_string(cluster),
                opacity=0.5
            ),
            row=row_hist, col=2
        )
    fig.update_xaxes(title_text="Decay Rate (slope)", row=row_hist, col=2)
    fig.update_yaxes(title_text="Count", row=row_hist, col=2)

    # Calculate stats for growth rates
    table_data_growth = []
    for cluster in clusters:
        mean = np.mean(slopes_mcmc[cluster], axis=0)
        lower_ci = np.percentile(slopes_mcmc[cluster], 2.5, axis=0)
        upper_ci = np.percentile(slopes_mcmc[cluster], 97.5, axis=0)
        table_data_growth.append([f'Cluster {cluster}', f'{mean:.4f}', f'{lower_ci:.4f} to {upper_ci:.4f}'])

    # Calculate stats for decay rates
    table_data_decay = []
    for cluster in clusters:
        mean_decay = np.mean(slopes_mcmc_decay[cluster], axis=0)
        lower_ci_decay = np.percentile(slopes_mcmc_decay[cluster], 2.5, axis=0)
        upper_ci_decay = np.percentile(slopes_mcmc_decay[cluster], 97.5, axis=0)
        table_data_decay.append(
            [f'Cluster {cluster}', f'{mean_decay:.4f}', f'{lower_ci_decay:.4f} to {upper_ci_decay:.4f}'])

    # Add growth rate table (second last row)
    fig.add_trace(
        go.Table(
            header=dict(
                values=['Cluster', 'Mean Growth Rate', '95% CI'],
                fill_color='lightblue',
                align='left'
            ),
            cells=dict(
                values=list(zip(*table_data_growth)),
                fill_color='white',
                align='left'
            )
        ),
        row=num_rows - 1, col=1
    )

    # Add decay rate table (last row)
    fig.add_trace(
        go.Table(
            header=dict(
                values=['Cluster', 'Mean Decay Rate', '95% CI'],
                fill_color='lightgreen',
                align='left'
            ),
            cells=dict(
                values=list(zip(*table_data_decay)),
                fill_color='white',
                align='left'
            )
        ),
        row=num_rows, col=1)

    # Update layout
    fig.update_layout(
        title="CLL Subclone Analysis",
        showlegend=True,
        height=400 * num_rows,  # Adjust height based on the number of rows
    )

    #     print(slopes_mcmc_decay)
    # Return HTML content
    return fig.to_html(full_html=False)


