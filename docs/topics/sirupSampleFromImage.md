# sirupSampleFromImage

## Get the sample name from an image label

## Description

`sirupSampleFromImage` returns the sample name in an image label.

## Usage

``` R
sirupSampleFromImage(image)
```

## Arguments

`image`

:   Character vector of image labels.

## Details

A multi-sample SIRUP object labels each spatial image as `"<sample>::<image>"`, so that image names stay unique across sections. The function returns the sample part of such a label. For a label without a sample part, it returns `NA`.

## Value

A character vector of sample names, with `NA` for each label that has no sample part.

## Examples

``` R
sirupSampleFromImage(c("sample1::image1", "image2"))
```
