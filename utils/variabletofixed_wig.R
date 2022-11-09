library(rtracklayer)

wig = import("/home/dhthutrang/test/test.wig")
print(wig)
export(wig, "/home/dhthutrang/test/test.wig_2", format="wig", dataFormat="fixedStep")