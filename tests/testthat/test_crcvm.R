
library(sf)


set.seed(123)


# Create a mock dataset for testing
mock_data <- data.frame(
  time = rep(1:10, 3),
  ID = factor(rep(1:3, each=10)))


# Loop to generate separate trajectories for each ID
for (i in 1:3) {
    mock_data[mock_data$ID == i, "x"] <- cumsum(rnorm(10))  
    mock_data[mock_data$ID == i, "y"] <- cumsum(rnorm(10))  
}

# Define a mock border (a simple square polygon for example)
large_rectangle_coords <- matrix(c(-10, -10, -10, 10, 10, 10, 10, -10,-10,-10), ncol = 2, byrow = TRUE)
small_rectangle_coords <- matrix(c(-8, -8, -8, 8, 8, 8, 8, -8,-8,-8), ncol = 2, byrow = TRUE)

# Create an sf object representing the rectangles
large_rectangle <- st_sfc(st_polygon(list(large_rectangle_coords)))
large_rectangle<-st_sf(geometry = large_rectangle)

small_rectangle <- st_sfc(st_polygon(list(small_rectangle_coords)))
small_rectangle<-st_sf(geometry = small_rectangle)

#change its crs
border<- st_sym_difference(large_rectangle, small_rectangle)



result <- get_BoundaryMetrics(mock_data, response=c("x","y"), border)


n_step=10

mock_data<-cbind(mock_data,result)

interpolation_data<-interpolate_BoundaryMetrics(mock_data,response=c("x","y"),border,n_step=n_step)

test_that("Interpolation has the right number of time steps",{
    
    expect_true(ncol(interpolation_data$BoundaryDistance)==n_step)
    expect_true(ncol(interpolation_data$BoundaryAngle)==n_step)
    
})


formulas=list("tau"=~1,"nu"=~1,"a"=~1,"b"=~1,"Dr"=~1,"Da"=~1,"sigma_D"=~1,
              "sigma_theta"=~1)

par0<-rep(1,8)


crcvm_interpolation_test1<-SDE$new(formulas = formulas,data=mock_data, type="CRCVM_SSM", response=c("x","y"),
                                   par0 = par0)


crcvm_interpolation_test1$setup()



BoundaryDistance_matrix=matrix(rep(mock_data$BoundaryDistance,each=10),ncol=10,byrow=TRUE)
BoundaryAngle_matrix=matrix(rep(mock_data$BoundaryAngle,each=10),ncol=10,byrow=TRUE)
crcvm_interpolation_test2<-SDE$new(formulas = formulas,data=mock_data, type="CRCVM_SSM", response=c("x","y"),
                             par0 = par0,other_data=list("interpolated_distance"=BoundaryDistance_matrix,
                                                         "interpolated_angle"=BoundaryAngle_matrix))

crcvm_interpolation_test2$setup()

crcvm_interpolation<-SDE$new(formulas = formulas,data=mock_data, type="CRCVM_SSM", response=c("x","y"),
                             par0 = par0,other_data=list("interpolated_distance"=interpolation_data$BoundaryDistance,
                                                          "interpolated_angle"=interpolation_data$BoundaryAngle))


crcvm_interpolation$setup()

test_that("Covariance and mean matrices are computed correctly", {
    


        expect_equal(round(crcvm_interpolation_test1$tmb_obj()$report()$T_matrices[,,7],7),
        round(crcvm_interpolation_test2$tmb_obj()$report()$T_matrices[,,7],7))


        expect_equal(round(crcvm_interpolation_test1$tmb_obj()$report()$Q_matrices[,,7],7),
        round(crcvm_interpolation_test2$tmb_obj()$report()$Q_matrices[,,7],7))
})



crcvm_interpolation_test1$fit()
crcvm_interpolation_test2$fit()
crcvm_interpolation$fit()
    

test_that("Fitting of CRCVM works correctly", {
    
    expect_equal(crcvm_interpolation_test1$tmb_obj()$report()$nllk,
                 crcvm_interpolation_test2$tmb_obj()$report()$nllk)
    
})

