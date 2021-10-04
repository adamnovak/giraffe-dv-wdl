# Giraffe/DeepVariant WDL Scripts

These WDL scripts are meant for testing running the Giraffe short read to graph
mapper and the DeepVariant variant caller together, in order to build an
accurate variant calling pipeline from the two.

`giraffe_and_deepvariant.wdl` runs the pipeline for a sample and evaluates the
resulting VCF against a truth set.

`generate_training_data.wdl` runs Giraffe and the indel realignment against a
series of samples, to generate DeepVariant training data.
