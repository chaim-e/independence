# Copyright (c) 2020 Chaim Even-Zohar
#
# This file is part of the R package "independence".
#
# "independence" is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# "independence" is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with "independence".  If not, see <https://www.gnu.org/licenses/>.


library(independence)

test_that("relative.order works",
{
    # input types
    expect_error(relative.order())
    expect_error(relative.order(1:100))
    expect_error(relative.order("abc", "def"))
    expect_error(relative.order(1,2,3,4))

    # length match
    expect_error(relative.order(1:100,1:101))
    expect_error(relative.order(1:101,1:100))

    # na.rm
    expect_error(relative.order(c(1,2,3,NA,NA,NA),c(NA,NA,NA,1,2,3)))
    expect_silent(relative.order(c(1,2,3,NA,NA,NA),c(NA,NA,NA,1,2,3),F,F))
    expect_length(relative.order(c(1,2,3,NA,NaN,Inf,-Inf),c(1,NA,3,4,5,6,7),T),2)
    expect_length(relative.order(c(1,2,3,NA,NaN,Inf,-Inf),c(1,NA,3,4,5,6,7),F),7)

    # collisions
    expect_warning(relative.order(c(1,1,1),c(2,2,3),T,T))
    expect_silent(relative.order(c(1,1,1),c(2,2,3),T,F))
    expect_warning(relative.order(as.integer(runif(10,0,9)),1:10))

    # function's function
    expect_equal(relative.order(1:10,1:10),1:10)
    expect_equal(relative.order(4:6,9:7),3:1)

    set.seed(1)
    expect_equal(relative.order(rnorm(7),rnorm(7)),c(3,2,6,5,4,1,7))
    x <- rnorm(100)
    y <- rnorm(100)
    expect_equal(relative.order(x, 1:100), order(x))
    expect_equal(relative.order(1:100, x), order(order(x)))
    expect_equal(relative.order(x,y), order(x[order(y)]))
    expect_equal(relative.order(x,y), relative.order(order(order(x)),y))
    expect_equal(relative.order(x,y), order(relative.order(y,x)))
})

test_that("rcpp functions work",
{
    # return types
    expect_is(independence:::max_hoeffding(), "numeric")
    expect_is(independence:::max_taustar(), "numeric")
    expect_is(independence:::calc_taustar(0:9), "numeric")
    expect_is(independence:::calc_hoeffding(0:9), "numeric")
    expect_is(independence:::calc_refined(0:9), "numeric")
    expect_type(independence:::calc_taustar(0:9), "double")
    expect_type(independence:::calc_hoeffding(0:9), "double")
    expect_type(independence:::calc_refined(0:9), "double")

    # selected values
    expect_equal(independence:::calc_taustar(c(2,1,0)),-1)
    expect_equal(independence:::calc_taustar(c(0,1,2,3)),2/3)
    expect_equal(independence:::calc_taustar(c(1,0,2,3)),2/3)
    expect_equal(independence:::calc_taustar(c(1,3,2,0)),-1/3)
    expect_equal(independence:::calc_taustar(c(1,2,0,3)),-1/3)
    expect_equal(independence:::calc_taustar(c(0,2,4,1,3,5)),0)
    expect_equal(independence:::calc_taustar(c(9,6,8,0,4,2,7,3,1,5)),1/42)
    expect_equal(independence:::calc_taustar(0:99),2/3)
    expect_equal(independence:::calc_taustar((0:100*24)%%101),-10/1111)
    expect_equal(independence:::calc_hoeffding(c(3,2,1,0)),-1)
    expect_equal(independence:::calc_hoeffding(c(0,1,2,4,3)),1/30)
    expect_equal(independence:::calc_hoeffding(c(1,4,2,3,0)),-1/60)
    expect_equal(independence:::calc_hoeffding(c(1,2,3,4,0)),0)
    expect_equal(independence:::calc_hoeffding(c(1,2,0,3,4)),0)
    expect_equal(independence:::calc_hoeffding(c(0,5,2,3,1,4)),-1/90)
    expect_equal(independence:::calc_hoeffding(c(9,6,8,0,4,2,7,3,1,5)),-1/945)
    expect_equal(independence:::calc_hoeffding(0:99),1/30)
    expect_equal(independence:::calc_refined(c(3,2,1,0)),-1)
    expect_equal(independence:::calc_refined(c(0,1,2,4,3)),1/90)
    expect_equal(independence:::calc_refined(c(1,2,0,3,4)),1/90)
    expect_equal(independence:::calc_refined(c(1,4,2,3,0)),-1/180)
    expect_equal(independence:::calc_refined(c(1,2,3,4,0)),-1/180)
    expect_equal(independence:::calc_refined(c(1,3,0,5,4,2)),0)
    expect_equal(independence:::calc_refined(c(9,6,8,0,4,2,7,3,1,5)),23/15120)
    expect_equal(independence:::calc_refined(0:99),1/90)

    # random values
    set.seed(1)
    x <- order(rnorm(12345))-1
    Tn <- independence:::calc_taustar(x)
    Dn <- independence:::calc_hoeffding(x)
    Rn <- independence:::calc_refined(x)
    expect_equal(Tn, -0.0000081824294564, tolerance = 1e-10, scale = abs(Tn))
    expect_equal(Dn, -0.0000002278996886, tolerance = 1e-10, scale = abs(Dn))
    expect_equal(Rn, -0.0000002269847163, tolerance = 1e-10, scale = abs(Rn))
    expect_equal(Tn, 12*Dn+24*Rn, tolerance = 1e-15)
})

test_that("independence tests test",
{
    # valid input
    for (test in c(tau.star.test,
                   hoeffding.D.test,
                   hoeffding.refined.test))
    {
        expect_error(test(rnorm(22),rnorm(33)))
        expect_error(test(1:3,11:13))
        expect_error(test(1:10,1:10 * NA))
        expect_error(test(1:10,1:10 * NaN))
        expect_error(test(1:10,1:10 * Inf))
        expect_warning(test(c(1,2,3,4,5,6),c(5,6,7,8,9,9)))
        expect_silent(test(1:6,1:6))
        if (requireNamespace("TauStar"))
        {
            expect_equal(is.na(test(1:6,1:6,precision = 0.0)$p.value), T)
            expect_equal(is.na(test(1:6,1:6,precision = 0.1)$p.value), F)
            expect_equal(is.na(test(1:6,1:6,precision = 1.0)$p.value), T)
        }
    }

    # independent case
    set.seed(2)
    xs <- rnorm(12345)
    ys <- rnorm(12345)
    ts <- tau.star.test(xs, ys)
    hd <- hoeffding.D.test(xs, ys)
    hr <- hoeffding.refined.test(xs, ys)
    for (it in list(ts,hd,hr))
    {
        expect_is(it, "indtest")
        expect_equal(it$n, 12345)
        if (requireNamespace("TauStar"))
        {
            expect_equal(it$p.value,
                         TauStar::pHoeffInd(it$scaled, lower.tail = F)[1])
            expect_gt(it$p.value, 0.001)
        }
    }
    expect_equal(ts$scaled, ts$Tn * 12344)
    expect_equal(hd$scaled, hd$Dn * 12344 * 36)
    expect_equal(hr$scaled, hr$Rn * 12344 * 36)

    # dependent case
    set.seed(3)
    xs <- rnorm(12345,0,1:12345+3000)
    ys <- rnorm(12345,0,1:12345+4000)
    ts <- tau.star.test(xs, ys)
    hd <- hoeffding.D.test(xs, ys)
    hr <- hoeffding.refined.test(xs, ys)
    for (it in list(ts,hd,hr))
    {
        expect_is(it, "indtest")
        expect_equal(it$n, 12345)
        if (requireNamespace("TauStar"))
        {
            expect_equal(it$p.value,
                         TauStar::pHoeffInd(it$scaled, lower.tail = F)[1])
            expect_lt(it$p.value, 0.001)
        }
    }
})
