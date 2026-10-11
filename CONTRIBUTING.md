# Contributing

Thanks for helping improve scAtlasLoom. Use [Issues](https://github.com/diediray/scatlasloom/issues/new/choose) for reproducible bugs, documentation fixes, and proposals for new figure modules.

## Report a rendering problem

Include the command, R and package versions (`sessionInfo()`), expected behavior, and the actual result. Prefer a small synthetic job with the same dimensions and column names. Remove patient identifiers, credentials, and private file paths before posting.

## Propose a change

Keep changes focused. Explain the research or usability need, identify the affected input contract, and preserve existing scientific review requirements. New modules should define their statistical unit, value scale, input alignment, and missing-data behavior. A new palette or layout should not change cell memberships or numerical results.

For R rendering changes, run the synthetic example and inspect PNG, PDF, and SVG outputs. Include a preview and report the platform used. New analysis methods need their own scientific rationale and relevant validation; a successful render alone is insufficient.
