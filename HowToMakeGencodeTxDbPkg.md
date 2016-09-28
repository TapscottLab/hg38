# How to Make a TxDb for Genome Build hg38
The aim of this document is to show how to make Bioconductor Transcripts database (TxDb) and package for Gencode track. Generally there are many ways to do it, and here I like to show how to extrack the current Gencode track from UCSC genome browser, which might be tricky, especially when Bionconductor 

## Build Gencode TxDb From UCSC Genome Brower
Let's start with loading two major libraries and open a browser session on R workspace.


```r
library(rtracklayer)
library(GenomicFeatures)
session <- browserSession()
genome(session) <- "hg38"
```

To check the track names and table names of Gencode on UCSC Genome Browser, you can do

```r
rtracklayer::trackNames(session)
```

Currently on on UCSC Genome Brower (9/26/2016), the track name  for Gencode is _Gencdoe v24_ and the table name is _knownGene_.  If you have tried


```r
txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene")
```
and if it works. Then you can ignore the rest.

**If it doesn't work**, that means there are some errors inhibit in the *GenomicFeatures* pakcage. You will need to take a detour and can use the following code as an example.


```r
tablename <- "knownGene"
track <- "GENCODE v24"
ucsc_txtable <-
    getTable(ucscTableQuery(session, track=track, table=tablename))
mapdef <- GenomicFeatures:::.howToGetTxName2GeneIdMapping(tablename)
txname2geneid <- GenomicFeatures:::.fetchTxName2GeneIdMappingFromUCSC(session, 
           track, tablename, mapdef)
```

```
## Download the knownToLocusLink table ...
```

```
## OK
```

```r
transcript_ids = NULL
txdb <- GenomicFeatures:::.makeTxDbFromUCSCTxTable(ucsc_txtable,
        txname2geneid$genes, 
        genome="hg38", tablename, track, txname2geneid$gene_id_type, 
        full_dataset = is.null(transcript_ids), circ_seqs = DEFAULT_CIRC_SEQS,
        taxonomyId = NA, miRBaseBuild = NA)
```

```
## Extract the 'transcripts' data frame ... OK
```

```
## Extract the 'splicings' data frame ...
```

```
## Warning in .extractCdsLocsFromUCSCTxTable(ucsc_txtable, exon_locs): UCSC data anomaly in 20500 transcript(s): the cds cumulative
##   length is not a multiple of 3 for transcripts 'uc057axw.1'
##   'uc031pko.2' 'uc057ayh.1' 'uc057ayr.1' 'uc057azl.1' 'uc057azn.1'
##   'uc057azp.1' 'uc057azt.1' 'uc057azy.1' 'uc057bad.1' 'uc057bap.1'
##   'uc057bav.1' 'uc057bay.1' 'uc057bbh.1' 'uc057bbv.1' 'uc057bcm.1'
##   'uc057bcp.1' 'uc057bcw.1' 'uc057bcx.1' 'uc057bdf.1' 'uc057bdj.1'
##   'uc057bdk.1' 'uc057bdl.1' 'uc057bdp.1' 'uc057beb.1' 'uc057beq.1'
##   'uc057bfn.1' 'uc057bfs.1' 'uc057bgm.1' 'uc057bgo.1' 'uc057bgt.1'
##   'uc057bhd.1' 'uc057bhi.1' 'uc057bio.1' 'uc057biu.1' 'uc057bji.1'
##   'uc057bjj.1' 'uc057bjk.1' 'uc057bkd.1' 'uc057bkf.1' 'uc057bkj.1'
##   'uc057bkl.1' 'uc057bkt.1' 'uc057bld.1' 'uc057bli.1' 'uc057blj.1'
##   'uc057blz.1' 'uc057bmf.1' 'uc057bmr.1' 'uc057bmv.1' 'uc057bni.1'
##   'uc057bnj.1' 'uc057bnk.1' 'uc057bnl.1' 'uc057bnv.1' 'uc057bob.1'
##   'uc057bod.1' 'uc057boe.1' 'uc057bpl.1' 'uc057bpm.1' 'uc057bpr.1'
##   'uc057bpw.1' 'uc057bqa.1' 'uc057bqb.1' 'uc057bqt.1' 'uc057bqx.1'
##   'uc057brc.1' 'uc057brd.1' 'uc057brp.1' 'uc057bsk.1' 'uc057bsq.1'
##   'uc057bsr.1' 'uc057bss.1' 'uc057bsv.1' 'uc057bsy.1' 'uc057bsz.1'
##   'uc057btc.1' 'uc057btg.1' 'uc057bti.1' 'uc057btj.1' 'uc057bua.1'
##   'uc057buh.1' 'uc284mmn.1' 'uc057buk.2' 'uc057bun.1' 'uc057bvb.1'
##   'uc057bvn.1' 'uc057bvo.1' 'uc057bvp.1' 'uc057bvq.1' 'uc057bvy.1'
##   'uc057bwe.1' 'uc057bwj.1' 'uc057bwu.1' 'uc057bwv.1' 'uc057bxk.1'
##   'uc057bxl.1' 'uc057byc.1' 'uc057byz.1' 'uc057bzv.1' 'uc057bzy.1'
##   'uc057cad.1' 'uc057cae.1' 'uc057caf.1' 'uc057cag.1' 'uc057cau.1'
##   'uc057caw.1' 'uc057cbs.1' 'uc057cby.1' 'uc057cca.1' 'uc284mmz.1'
##   'uc057cco.1' 'uc057ccq.1' 'uc057cdi.1' 'uc057cdn.1' 'uc057cdo.1'
##   'uc057cdp.1' 'uc057cds.1' 'uc057cdw.1' 'uc057cdz.1' 'uc057cee.1'
##   'uc057cef.1' 'uc057cei.1' 'uc057cej.1' 'uc057cel.1' 'uc057cem.1'
##   'uc057cen.1' 'uc057ceo.1' 'uc057cep.1' 'uc057ceq.1' 'uc057cfd.1'
##   'uc057cfg.1' 'uc057cgj.1' 'uc057cgl.1' 'uc057cgn.1' 'uc057cgq.1'
##   'uc057cgt.1' 'uc057cgu.1' 'uc057chi.1' 'uc057chk.1' 'uc057chl.1'
##   'uc057chq.1' 'uc057chr.1' 'uc057chs.1' 'uc057cht.1' 'uc057chu.1'
##   'uc057cil.1' 'uc057cjn.1' 'uc001aue.2' 'uc057clh.1' 'uc057clq.1'
##   'uc057clz.1' 'uc057cmf.1' 'uc057cmg.1' 'uc057cmw.1' 'uc057cna.1'
##   'uc057cnf.1' 'uc057cnm.1' 'uc057cno.1' 'uc057cnp.1' 'uc057cnr.1'
##   'uc057cns.1' 'uc057cnu.1' 'uc057cny.1' 'uc057coi.1' 'uc057cox.1'
##   'uc057coz.1' 'uc057cpa.1' 'uc057cpc.1' 'uc057cpd.1' 'uc057cpj.1'
##   'uc057cpn.1' 'uc057cpo.1' 'uc057cqx.1' 'uc057cqy.1' 'uc057crr.1'
##   'uc057cso.1' 'uc057csp.1' 'uc057csr.1' 'uc057ctm.1' 'uc057ctn.1'
##   'uc057cto.1' 'uc057cva.1' 'uc057cvc.1' 'uc057cve.1' 'uc057cvo.1'
##   'uc057cvp.1' 'uc057cvx.1' 'uc057cwh.1' 'uc057cwm.1' 'uc057cxh.1'
##   'uc057cxk.1' 'uc057cxl.1' 'uc057cxm.1' 'uc057cxn.1' 'uc057cxo.1'
##   'uc057cyx.1' 'uc057cyy.1' 'uc057czk.1' 'uc057czs.1' 'uc057czt.1'
##   'uc057czv.1' 'uc057dad.1' 'uc057dah.1' 'uc057dai.1' 'uc057dap.1'
##   'uc057daz.1' 'uc057dbp.1' 'uc057dbx.1' 'uc057dby.1' 'uc057dce.1'
##   'uc057dcf.1' 'uc057dcg.1' 'uc057dcj.1' 'uc057dcr.2' 'uc057dcu.1'
##   'uc057dda.1' 'uc057ddn.1' 'uc057ddq.1' 'uc057ddt.1' 'uc057ddw.1'
##   'uc057ddz.1' 'uc057deb.1' 'uc057dek.1' 'uc057deo.1' 'uc057dep.1'
##   'uc057des.1' 'uc057dfu.1' 'uc057dfy.1' 'uc057dgr.1' 'uc057dgx.1'
##   'uc057dhe.1' 'uc057dhh.1' 'uc057dhp.1' 'uc057dhw.1' 'uc057diw.1'
##   'uc057djc.1' 'uc057djl.1' 'uc057djo.1' 'uc057djs.1' 'uc057dkb.1'
##   'uc057dlb.1' 'uc057dlc.1' 'uc057dln.1' 'uc057dme.1' 'uc057dnf.1'
##   'uc057dng.1' 'uc057dni.1' 'uc057dnj.1' 'uc057dnr.1' 'uc057dnt.1'
##   'uc057dof.1' 'uc057dol.1' 'uc057dom.1' 'uc057dpl.1' 'uc057dps.1'
##   'uc057dqf.1' 'uc057dqg.1' 'uc057dqi.1' 'uc057dqj.1' 'uc057dqk.1'
##   'uc057dqr.1' 'uc057dqw.1' 'uc057dqx.1' 'uc057dqy.1' 'uc057dqz.1'
##   'uc057drb.1' 'uc057drj.1' 'uc057drk.1' 'uc057drl.1' 'uc057drq.1'
##   'uc057drx.1' 'uc057dsn.1' 'uc057dtd.1' 'uc057dtg.1' 'uc057dth.1'
##   'uc057dti.1' 'uc057dtl.1' 'uc057dtm.1' 'uc057dtp.1' 'uc057dtw.1'
##   'uc057dtx.1' 'uc057duc.1' 'uc057duh.1' 'uc057dva.1' 'uc057dvg.1'
##   'uc057dvh.1' 'uc057dvy.1' 'uc057dwt.1' 'uc057dwx.1' 'uc057dxb.1'
##   'uc057dxi.1' 'uc057dxt.1' 'uc284mnw.1' 'uc057dyk.1' 'uc057dym.1'
##   'uc057dyn.1' 'uc057dyo.1' 'uc057dzs.1' 'uc057ear.1' 'uc057ecs.1'
##   'uc057ect.1' 'uc057ecv.1' 'uc057ecw.1' 'uc057edb.1' 'uc057edj.1'
##   'uc057edk.1' 'uc057edy.1' 'uc057egg.1' 'uc031twa.2' 'uc057egt.1'
##   'uc057ehd.1' 'uc057ehe.1' 'uc057ehf.1' 'uc057eif.1' 'uc057eip.1'
##   'uc057ejn.1' 'uc057eju.1' 'uc057ejv.1' 'uc057ekb.1' 'uc057eks.1'
##   'uc057elc.1' 'uc057eld.1' 'uc057elf.1' 'uc057eli.1' 'uc057elz.1'
##   'uc057ema.1' 'uc057emb.1' 'uc057emd.1' 'uc057emg.1' 'uc057emh.1'
##   'uc057emt.1' 'uc057emw.1' 'uc057emx.1' 'uc057emz.1' 'uc057enb.1'
##   'uc057eow.1' 'uc057eqq.1' 'uc057eqz.1' 'uc057erc.1' 'uc057erk.1'
##   'uc057erv.1' 'uc057esb.1' 'uc057eso.1' 'uc057est.1' 'uc284moh.1'
##   'uc057etm.1' 'uc057ett.1' 'uc057eud.1' 'uc057euu.1' 'uc057eve.1'
##   'uc057evk.1' 'uc057evt.1' 'uc057ewd.1' 'uc057ewq.1' 'uc057ewy.1'
##   'uc057ewz.1' 'uc057exb.1' 'uc057exd.1' 'uc057exf.1' 'uc057exh.1'
##   'uc057exp.1' 'uc057ezt.1' 'uc057fah.1' 'uc057fam.1' 'uc057fao.1'
##   'uc057fas.1' 'uc057fau.1' 'uc057fav.1' 'uc057fbj.1' 'uc057fbm.1'
##   'uc057fbn.1' 'uc057fbq.1' 'uc057fbv.1' 'uc057fcd.1' 'uc057fcf.1'
##   'uc057fcp.1' 'uc057fdb.1' 'uc057fde.1' 'uc057fdj.1' 'uc057fdl.1'
##   'uc057fds.1' 'uc057fdt.1' 'uc057fdu.1' 'uc057fdw.1' 'uc057fdz.1'
##   'uc057fea.1' 'uc057fec.1' 'uc057fed.1' 'uc057fem.1' 'uc057fex.1'
##   'uc057fez.1' 'uc057fff.1' 'uc057ffp.1' 'uc057ffs.1' 'uc057fft.1'
##   'uc057fgl.1' 'uc057fgp.1' 'uc057fgq.1' 'uc057fgw.1' 'uc057fip.1'
##   'uc057fjo.1' 'uc057fjt.1' 'uc057fju.1' 'uc057fka.1' 'uc057fkf.1'
##   'uc057fla.1' 'uc057fld.1' 'uc057flr.1' 'uc057fmf.1' 'uc057fmh.1'
##   'uc057fnt.1' 'uc057fnv.1' 'uc057fnw.1' 'uc057fnx.1' 'uc057foh.1'
##   'uc057foj.1' 'uc057fok.1' 'uc057fol.1' 'uc057fpz.1' 'uc057fqb.1'
##   'uc057fqc.1' 'uc057fqd.1' 'uc057fqg.1' 'uc057frd.1' 'uc057frl.1'
##   'uc057fsd.1' 'uc057fsf.1' 'uc057fti.1' 'uc057ftk.1' 'uc057fty.1'
##   'uc057fua.1' 'uc057fug.1' 'uc057fux.1' 'uc057fuz.1' 'uc057fvc.1'
##   'uc057fve.1' 'uc057fvf.1' 'uc057fvi.1' 'uc057fvr.1' 'uc057fvx.1'
##   'uc057fww.1' 'uc057fwy.1' 'uc057fxp.1' 'uc057fxy.1' 'uc057fyh.1'
##   'uc057fyo.1' 'uc057fzw.1' 'uc057fzx.1' 'uc057fzz.1' 'uc057gac.1'
##   'uc057gar.1' 'uc057gas.1' 'uc057gau.1' 'uc057gaw.1' 'uc057gbb.1'
##   'uc057gbd.1' 'uc057gcn.1' 'uc057gdi.1' 'uc057gdj.1' 'uc057gdk.1'
##   'uc057gdm.1' 'uc057gdn.1' 'uc057gdq.1' 'uc057gdr.1' 'uc057geh.1'
##   'uc057gel.1' 'uc057gev.1' 'uc057gfe.1' 'uc057gfk.1' 'uc057gfz.1'
##   'uc057ggb.1' 'uc057gge.1' 'uc057ggh.1' 'uc057ggr.1' 'uc057ghr.1'
##   'uc057ghv.1' 'uc057gio.1' 'uc057gis.1' 'uc057giw.1' 'uc057gjg.1'
##   'uc057gks.1' 'uc057gkv.1' 'uc057gkx.1' 'uc057gky.1' 'uc057gla.1'
##   'uc057gll.1' 'uc057glr.1' 'uc057glw.1' 'uc057gmb.1' 'uc057gmh.1'
##   'uc057gml.1' 'uc057gnx.1' 'uc057goa.1' 'uc057goh.1' 'uc057gok.1'
##   'uc057gol.1' 'uc057gom.1' 'uc057gon.1' 'uc057gos.1' 'uc057gqk.1'
##   'uc057gql.1' 'uc057gqm.1' 'uc057gri.1' 'uc057grn.1' 'uc057gro.1'
##   'uc057gsy.1' 'uc057gsz.1' 'uc057gtf.1' 'uc057gtg.1' 'uc057gth.1'
##   'uc057gub.1' 'uc057gvw.1' 'uc057gyq.1' 'uc057gys.1' 'uc057gyu.1'
##   'uc057gyv.1' 'uc057gyw.1' 'uc057gyx.1' 'uc057hab.1' 'uc057hac.1'
##   'uc057has.1' 'uc057hbk.1' 'uc057hbp.1' 'uc057hby.1' 'uc057hcd.1'
##   'uc057hcg.1' 'uc057hci.1' 'uc284mpr.1' 'uc284mpv.1' 'uc057hda.1'
##   'uc284mqe.1' 'uc284mqg.1' 'uc057hdl.1' 'uc057hdm.1' 'uc057hdn.1'
##   'uc057hfu.1' 'uc057hfy.1' 'uc057hhq.1' 'uc057hhs.1' 'uc057hio.1'
##   'uc057his.1' 'uc057hjh.1' 'uc057hkl.1' 'uc057hkm.1' 'uc057hkn.1'
##   'uc057hko.1' 'uc057hkq.1' 'uc057hks.1' 'uc057hlb.1' 'uc057hlc.1'
##   'uc057hld.1' 'uc057hle.1' 'uc057hpk.1' 'uc057hpu.1' 'uc057hpy.1'
##   'uc057hqc.1' 'uc057hqd.1' 'uc057hqe.1' 'uc057hqs.1' 'uc057hqt.1'
##   'uc057hqv.1' 'uc057htb.1' 'uc057htj.1' 'uc057htt.1' 'uc057htu.1'
##   'uc057htz.1' 'uc057hud.1' 'uc057huj.1' 'uc057hum.1' 'uc057hun.1'
##   'uc057hup.1' 'uc057hur.1' 'uc057hus.1' 'uc057hut.1' 'uc057hvu.1'
##   'uc057hwk.1' 'uc057hxj.1' 'uc057hxl.1' 'uc057hyq.1' 'uc057hzv.1'
##   'uc057iau.1' 'uc057iaw.1' 'uc057iaz.1' 'uc057ibb.1' 'uc057ibk.1'
##   'uc057icd.2' 'uc057ice.2' 'uc057icw.1' 'uc057icy.1' 'uc057idi.1'
##   'uc057ido.1' 'uc057iem.1' 'uc057ien.1' 'uc057ifb.1' 'uc057ifc.1'
##   'uc057ifd.1' 'uc057ife.1' 'uc057iff.1' 'uc057ify.1' 'uc057igc.1'
##   'uc057igd.1' 'uc057igf.1'
```

```
## OK
```

```
## Download and preprocess the 'chrominfo' data frame ...
```

```
## OK
```

```
## Prepare the 'metadata' data frame ...
```

```
## OK
```

```
## Make the TxDb object ...
```

```
## OK
```

```r
txdb
```

```
## TxDb object:
## # Db type: TxDb
## # Supporting package: GenomicFeatures
## # Data source: UCSC
## # Genome: hg38
## # Organism: Homo sapiens
## # Taxonomy ID: 9606
## # UCSC Table: knownGene
## # UCSC Track: GENCODE v24
## # Resource URL: http://genome.ucsc.edu/
## # Type of Gene ID: Entrez Gene ID
## # Full dataset: yes
## # miRBase build ID: NA
## # transcript_nrow: 197782
## # exon_nrow: 581036
## # cds_nrow: 293052
## # Db created by: GenomicFeatures package from Bioconductor
## # Creation time: 2016-09-28 11:51:21 -0700 (Wed, 28 Sep 2016)
## # GenomicFeatures version at creation time: 1.24.5
## # RSQLite version at creation time: 1.0.0
## # DBSCHEMAVERSION: 1.1
```

Once you have a TxDb instance, you can build a TxDb package. Below is an example:

```r
makeTxDbPackage(txdb, version = "1.24.0",
                maintainer = "Chao-Jen Wong <cwon2@fredhutch.org>",
                author = "Chao-Jen Wong", destDir = "~/tapscott/hg38")
```
				
## Ensembl From BioMart
Want to use the most recent annotation on Ensembl? See example below:

```r
library(biomaRt)
library(GenomicFeatures)
listMarts(host="www.ensembl.org")
datasets <- listDatasets(biomaRt::useMart(
    biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
ensembl.txdb <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                                    dataset="hsapiens_gene_ensembl",
                                    host="www.ensembl.org")
makeTxDbPackage(ensembl.txdb, version="1.85.0", author="Chao-Jen Wong",
                maintainer="Chao-Jen Wong <cwon2@fredhutch.org>")
```
				
## Session Information

```r
sessionInfo()
```

```
## R Under development (unstable) (2016-09-25 r71358)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.3 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] GenomicFeatures_1.24.5 AnnotationDbi_1.34.4   Biobase_2.32.0        
##  [4] rtracklayer_1.32.2     GenomicRanges_1.24.2   GenomeInfoDb_1.8.7    
##  [7] IRanges_2.6.1          S4Vectors_0.10.3       BiocGenerics_0.18.0   
## [10] knitr_1.14             BiocInstaller_1.22.3  
## 
## loaded via a namespace (and not attached):
##  [1] XVector_0.12.1             magrittr_1.5              
##  [3] zlibbioc_1.18.0            GenomicAlignments_1.8.4   
##  [5] BiocParallel_1.6.6         stringr_1.1.0             
##  [7] tools_3.4.0                SummarizedExperiment_1.2.3
##  [9] DBI_0.5-1                  formatR_1.4               
## [11] bitops_1.0-6               biomaRt_2.28.0            
## [13] RCurl_1.95-4.8             RSQLite_1.0.0             
## [15] evaluate_0.9               stringi_1.1.1             
## [17] Biostrings_2.40.2          Rsamtools_1.24.0          
## [19] XML_3.98-1.4
```
