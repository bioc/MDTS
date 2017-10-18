context('Test if all necessary aspects of the package was installed correctly')

test_that('Testing if MDTS is installed and loading properly', {
	expect_silent(library(MDTS))
})

base_path = system.file('extdata', package='MDTS')
setwd(base_path)
test_that('Testing if sample data can be loaded', {
	expect_silent(load('bins.RData'))
	expect_silent(load('counts.RData'))
	expect_silent(load('pD.RData'))
})




