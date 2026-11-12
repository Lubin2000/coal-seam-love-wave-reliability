param(
    [switch]$SkipWavefields,
    [switch]$FiguresOnly
)

$ErrorActionPreference = "Stop"
$ExperimentRoot = Join-Path $PSScriptRoot "analysis"

$EligibleExperiments = @(
    "e13_acquisition_information.py",
    "e14_mode_availability.py",
    "e15_automatic_mode_picking.py",
    "e16_robust_mode_tracking.py",
    "e17_robust_likelihood.py",
    "e18_robust_posterior_calibration.py",
    "e25_qedispinv_validation.py",
    "e26_additional_picker_baseline.py",
    "e27_bootstrap_uncertainty.py"
)

$WavefieldExperiments = @(
    "e22_love_forward_validation.py",
    "e24_wavefield_tracking_inversion.py",
    "e28_sls_attenuation_validation.py"
)

if (-not $FiguresOnly) {
    $Experiments = $EligibleExperiments
    if (-not $SkipWavefields) {
        $Experiments += $WavefieldExperiments
    }
    foreach ($Script in $Experiments) {
        Write-Host "Running $Script"
        & py (Join-Path $ExperimentRoot $Script)
        if ($LASTEXITCODE -ne 0) {
            throw "$Script failed with exit code $LASTEXITCODE"
        }
    }
}

& py (Join-Path $PSScriptRoot "paper\make_figure6_detection.py")
if ($LASTEXITCODE -ne 0) {
    throw "Figure 6 detection generation failed with exit code $LASTEXITCODE"
}

& py (Join-Path $PSScriptRoot "paper\make_figures.py")
if ($LASTEXITCODE -ne 0) {
    throw "Figure generation failed with exit code $LASTEXITCODE"
}

Write-Host "Core paper results and figures are up to date."
