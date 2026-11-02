# Resume an interrupted orchestration

Continue the existing migration. Do not restart it, recreate repositories or reinterpret previous completed work as permission for broader actions.

1. Read START_HERE, shared rules, stack.dcf, the evidence STATUS.dcf, phase-status.csv, latest handoffs, frozen contract and current user instructions.
2. Reconcile local paths, branches, worktree states and actual HEADs with the ledger. Preserve dirty work. A recorded phase is not proof its commit is still the active candidate.
3. Inspect any in-flight worker or tool state. Do not assign a second writer to a file set already in use. Resume or take over only after ownership is clear.
4. For an uncertain remote create/push/PR operation, inspect the server before retrying. Record actual owner, visibility, branches and SHAs. Avoid duplicate actions.
5. Determine the earliest incomplete or invalidated prerequisite in phases.csv. Reuse completed evidence when the relevant candidate code/contract is unchanged; rerun checks only for a concrete change or unmet gate.
6. If the contract changed, update dependent work and fixtures deliberately. Do not combine results from incompatible contract versions.
7. Continue all independent authorized work. Ask for user input only for an actual missing decision or access constraint that cannot be resolved from the record; identify the concrete action and reason.
8. Keep state current and resume the implementation/review loop through P12. The final report must distinguish completed local work, verified checks and actual remote status.

