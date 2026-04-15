test_that("story-local TC-9.4.19 replay summary keeps weakest-case seed offsets", {
  diagnostics <- lwdid:::.story_local_wcb_power_diagnostics(
    seed = 100L,
    replay_seeds = c(100L, 146L)
  )
  replay_summary <- diagnostics$result$targeted_replay_summary

  expect_type(replay_summary, "list")
  expect_equal(replay_summary$weakest_hit$seed, 100L, tolerance = 0)
  expect_equal(replay_summary$weakest_hit$attempt_id, 1L, tolerance = 0)
  expect_equal(replay_summary$weakest_hit$replay_seed_offset, 0L, tolerance = 0)
  expect_equal(replay_summary$weakest_miss$seed, 146L, tolerance = 0)
  expect_equal(replay_summary$weakest_miss$attempt_id, 47L, tolerance = 0)
  expect_equal(replay_summary$weakest_miss$replay_seed_offset, 46L, tolerance = 0)
})
