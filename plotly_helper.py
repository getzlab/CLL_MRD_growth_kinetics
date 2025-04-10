import plotly.graph_objects as go
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