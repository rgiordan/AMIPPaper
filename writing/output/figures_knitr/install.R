library(devtools)

inst_packages <- rownames(installed.packages())

ConditionalInstall <- function(package) {
    if (!(package %in% inst_packages)) {
        install.package(package)
    } else {
        print(sprintf("%s already installed.", package))
    }
}
ConditionalInstall("stargazer")
ConditionalInstall("AER")
ConditionalInstall("reticulate")

# Actually you may need to do additional work for torch...
ConditionalInstall("torch")

if (!("zaminfluence" %in% inst_packages)) {
    devtools::install_github("https://github.com/rgiordan/zaminfluence/",
                            ref="v0.4",
                            subdir="zaminfluence",
                            force=TRUE,
                            upgrade="never")
} else {
        print(sprintf("%s already installed.", "zaminfluence"))
}
