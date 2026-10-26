#!/usr/bin/env Rscript
# Copyright (c) 2026 Raban Heller. MIT licence; see ../../LICENSE.md.
# Public, installed-package inventory of the semantic claims and actual code.
args <- commandArgs(TRUE)
if(length(args)!=1L) stop("usage: Rscript tools/product/inventory-advisor.R OUTPUT_DIRECTORY")
out<-path.expand(args[1])
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE))) stop("output directory must be new or empty")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
library(inferguidR)
package<-find.package("inferguidR")
if(!file.exists(file.path(package,"Meta","package.rds"))) stop("installed package required")
expected<-Sys.getenv("INFERGUIDR_EXPECTED_LIBRARY")
if(nzchar(expected)&&normalizePath(file.path(expected,"inferguidR"))!=normalizePath(package)) stop("unexpected package library")
rs<-advisor_rules()
variants<-jsonlite::read_json(system.file("advisor","analysis-variants.json",package="inferguidR"))
crosswalk<-utils::read.csv(system.file("advisor","source-crosswalk.csv",package="inferguidR"),stringsAsFactors=FALSE,na.strings=character())
rows<-list(); detail<-list()
for(v in variants$variants) {
 id<-v$analysis_variant_id
 primary<-Filter(function(o) identical(o$analysis_variant_id,id),rs$outcomes)
 alternatives<-Filter(function(o) any(vapply(o$alternatives,function(a) identical(a$analysis_variant_id,id)&&isTRUE(a$override_allowed),TRUE)),rs$outcomes)
 rejected<-Filter(function(o) any(vapply(o$rejected_routes,function(a) identical(a$analysis_variant_id,id),TRUE)),rs$outcomes)
 path<-system.file("advisor","templates",paste0("tpl.",id,".v1.R"),package="inferguidR")
 text<-if(nzchar(path)) readLines(path,warn=FALSE,encoding="UTF-8") else character()
 calls<-sort(unique(unlist(regmatches(text,gregexpr("[A-Za-z][A-Za-z0-9.]*::[A-Za-z][A-Za-z0-9._]*",text)))))
 # Retain actual implementation declarations, not only the historical sketches.
 claims<-text[grepl("estimand|estimate|assumption|interpretation|limitation|method=|method =|alternative|conf.level|exact|continuity|symmetr|parallel|bootstrap",text,ignore.case=TRUE)]
 cw<-crosswalk[crosswalk$analysis_variant_id==id,,drop=FALSE]
 assumptions<-unlist(lapply(c(primary,alternatives),function(o) lapply(o$assumptions,function(a)
   list(outcome_id=o$outcome_id,assumption_id=a$assumption_id,text_en=a$text_en,text_de=a$text_de))),recursive=FALSE)
 available<-nzchar(path)
 status<-if(length(primary)) "primary_route" else if(length(alternatives)) "allowed_override_only" else "rejected_only"
 relation<-if(nrow(cw)) paste(sort(unique(cw$relationship)),collapse="; ") else "extension_no_direct_source_leaf"
 rows[[length(rows)+1L]]<-data.frame(variant=id,route=status,template=if(available) "available" else "absent",
   estimand=v$estimand_en,source_relation=relation,primary_outcomes=paste(names(primary),collapse="; "),
   allowed_override_outcomes=paste(names(alternatives),collapse="; "),
   assumptions=length(assumptions),actual_namespaced_calls=paste(calls,collapse="; "),
   review_scope="Executable teaching identity and documented assumptions; not general estimator validity or numerical certification")
 detail[[id]]<-list(variant_id=id,estimand=v$estimand_en,metadata_notes=v$notes,
   route=status,primary_outcomes=names(primary),allowed_override_outcomes=names(alternatives),
   rejected_outcomes=names(rejected),outcome_assumptions=assumptions,
   source_crosswalk=cw,source_scope=if(nrow(cw)) "Recorded comparison with bundled source-tree snapshots; source equality is not asserted" else
     "Package extension without a direct leaf in the two source trees; requires method-specific literature and review",
   template_sha256=if(available) digest::digest(file=path,algo="sha256") else NULL,
   implementation_declarations=claims,actual_namespaced_calls=calls,
   historical_call_sketches=v$r_functions,
   interpretation="Call sketches are historical metadata; actual template code and its assumptions govern the example. A parsed or tested script is not broad statistical validation.")
}
inventory<-do.call(rbind,rows)
utils::write.csv(inventory,file.path(out,"variants.csv"),row.names=FALSE)
jsonlite::write_json(detail,file.path(out,"variant-details.json"),auto_unbox=TRUE,pretty=TRUE,null="null",na="null")
summary<-list(variants=nrow(inventory),available_templates=sum(inventory$template=="available"),
 rejected_only=inventory$variant[inventory$route=="rejected_only"],
 source_extensions=inventory$variant[inventory$source_relation=="extension_no_direct_source_leaf"],
 rules_version=rs$versions$rules_version,package_path=package,package_version=as.character(utils::packageVersion("inferguidR")),
 source_sha256=Sys.getenv("INFERGUIDR_SOURCE_SHA256"),
 scope="Every bundled variant: actual template hash/declarations/calls, estimand, outcome assumptions and precise recorded source relationship. No claim of source parity or full semantic/numerical certification.")
jsonlite::write_json(summary,file.path(out,"summary.json"),auto_unbox=TRUE,pretty=TRUE,null="null")
if(nrow(inventory)!=69L || any(inventory$template=="absent" & inventory$route!="rejected_only")) stop("unexpected variant or executable-template inventory")
