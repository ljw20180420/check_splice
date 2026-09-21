// 1. Inject a styled scrollable container below the graph matrix
const graphDiv = document.getElementById('{{plot_id}}');
const container = document.createElement('div');
container.innerHTML = `
    <div style="margin-top: 20px; font-family: sans-serif; width: 100%; max-width: 800px;">
        <h4 style="margin-bottom: 5px;">Data Point Details</h4>
        <div id="scroll-window" style="height: 150px; overflow-y: auto; border: 1px solid #ccc; padding: 12px; background: #f9f9f9; border-radius: 4px; line-height: 1.5;">
            Click on a data point in the plot above to view its detailed text.
        </div>
    </div>
`;
graphDiv.after(container);

// 2. Load the descriptions data
const textArray = { js_descriptions };

// 3. Listen for the Plotly click event
graphDiv.on('plotly_click', function (data) {
    {
        if (data && data.points && data.points.length > 0) {
            {
                const pointIndex = data.points[0].pointIndex;
                const targetWindow = document.getElementById('scroll-window');

                // Update the scrollable box with the text
                targetWindow.innerHTML = textArray[pointIndex] || "No description available.";
                // Reset scroll position back to top on every click
                targetWindow.scrollTop = 0;
            }
        }
    }
});
