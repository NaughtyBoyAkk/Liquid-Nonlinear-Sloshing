function genGeom(k1::String,c1::String,k2::String,c2::String)
    # 读入自由面的几何信息
    fU=open(k1,"r")
    U1=readdlm(fU,Float64)
    U1=[U1[1]; U1[:]; U1[end]] # Rhino导出的knot点首尾各须增加一个重复度
    close(fU)
    #
    fc=open(c1,"r")
    cpts1=readdlm(fc,Float64)
    close(fc)
    w1 = cpts1[:,4]
    cpts1=cpts1[:,1:2]
    # 读入湿表面的几何信息
    fU=open(k2,"r")
    U2=readdlm(fU,Float64)
    U2=[U2[1]; U2[:]; U2[end]] # Rhino导出的knot点首尾各须增加一个重复度
    close(fU)
    #
    fc=open(c2,"r")
    cpts2=readdlm(fc,Float64)
    close(fc)
    w2 = cpts2[:,4]
    cpts2 = cpts2[:,1:2]
    #
    return U1,cpts1,w1,U2,cpts2,w2
end
