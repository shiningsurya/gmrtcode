# gmrtcode

| Name | Description and usecase |
|------|-------------------------|
| `gmrtfits` | Provides an interface to write search mode `PSRFITS` format files which has been made to work with uGMRT data formats |
| `singlepol8bit` | Filterbanks 8bit real baseband data of an antenna and stores the autocorrelation product in PSRFITS format |
| `dualpol8bit` | Filterbanks 8bit real baseband data of an antenna and stores the coherence products in PSRFITS format |
| `dualpol8bitgpu` | Does what `dualpol8bit` does, but FFTs and detection are done with GPUs |
| `foldgpu` | Detects, de-disperses and folds 8bit real baseband data |
| `foldvoltgpu` | Folds 8bit real baseband data |
