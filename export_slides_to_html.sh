#!/bin/bash

#jupyter nbconvert final_presentation_internship.ipynb --to slides --TemplateExporter.exclude_input=True --EmbedImagesPreprocessor.embed_images=True 

jupyter nbconvert final_presentation_internship.ipynb --to slides \
    --SlidesExporter.reveal_scroll=True \
    --no-input \
    --EmbedImagesPreprocessor.embed_images=True


# Check if the HTML file was created successfully
if [ -f "final_presentation_internship.slides.html" ]; then
    echo "HTML file created successfully!"

    # Open in the default web browser (Linux)
    xdg-open final_presentation_internship.slides.html
else
    echo "Failed to create HTML file."
fi



