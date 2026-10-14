#' inferguidR: Guided Inference, Power Analysis and Diagnostic Accuracy
#'
#' inferguidR plans and explains statistical analyses. Its core is a registry
#' of power-analysis procedures modelled on the test inventory of the
#' 'G*Power' 3.1 program, typed analysis specifications ([analysis_spec()]),
#' results that always carry a status ([calculate_power()]), canonical JSON
#' with SHA-256 hashes ([ig_to_json()], [ig_hash()]) and safe import
#' ([ig_from_json()]).
#'
#' @section Status of the implementation:
#' Engines are added method by method. A method without an engine returns a
#' result with status `"unsupported"` and no number. No agreement with
#' 'G*Power' is claimed: no application reference captures exist yet, so the
#' G*Power reference of every inventory method is `"reference_pending"`, and
#' extensions that are not G*Power procedures (Boschloo's test) have none
#' (`"not_applicable"`; see [validation_status()]).
#'
#' @section Options:
#' * `inferguidR.reference_archive_dir`: path of a local reference archive
#'   (development only; default unset).
#' * `inferguidR.validate_json_schema`: validate imported JSON against the
#'   bundled JSON Schemas when 'jsonvalidate' is installed (default `TRUE`).
#'
#' See [ig_config()] for the current values.
#'
#' @keywords internal
"_PACKAGE"
