!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exchange numerical matrix payloads and integer metadata between MPI ranks.
! **************************************************************************************************
MODULE matrix_communication
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'matrix_communication'
   PUBLIC :: matrix_buffer_type, allocate_matrix_buffers, exchange_matrix_buffers

   TYPE matrix_buffer_type
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: msg
      ! Metadata records are columns, contiguous within each rank's segment.
      INTEGER, DIMENSION(:, :), ALLOCATABLE              :: indx
      ! Zero-based segment offsets, including the final total at index nprocs.
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: msg_offset, indx_offset
   END TYPE matrix_buffer_type

CONTAINS

! **************************************************************************************************
!> \brief Allocate flat, rank-grouped matrix buffers and their zero-based offsets.
!> \param send_count numerical payload length sent to each rank
!> \param recv_count numerical payload length received from each rank
!> \param n_index number of integer coordinates per metadata record
!> \param buffer_send output send buffer
!> \param buffer_rec output receive buffer
!> \param send_blocks optional metadata-record counts, distinct from payload counts
!> \param recv_blocks optional received metadata-record counts
! **************************************************************************************************
   SUBROUTINE allocate_matrix_buffers(send_count, recv_count, n_index, buffer_send, buffer_rec, &
                                      send_blocks, recv_blocks)
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: send_count, recv_count
      INTEGER, INTENT(IN)                                :: n_index
      TYPE(matrix_buffer_type), INTENT(OUT)              :: buffer_send, buffer_rec
      INTEGER, DIMENSION(0:), INTENT(IN), OPTIONAL       :: send_blocks, recv_blocks

      INTEGER                                            :: nprocs, peer

      CPASSERT(PRESENT(send_blocks) .EQV. PRESENT(recv_blocks))
      nprocs = SIZE(send_count)
      ALLOCATE (buffer_send%msg_offset(0:nprocs), buffer_rec%msg_offset(0:nprocs), &
                buffer_send%indx_offset(0:nprocs), buffer_rec%indx_offset(0:nprocs))
      buffer_send%msg_offset(0) = 0
      buffer_rec%msg_offset(0) = 0
      buffer_send%indx_offset(0) = 0
      buffer_rec%indx_offset(0) = 0
      DO peer = 0, nprocs - 1
         buffer_send%msg_offset(peer + 1) = buffer_send%msg_offset(peer) + send_count(peer)
         buffer_rec%msg_offset(peer + 1) = buffer_rec%msg_offset(peer) + recv_count(peer)
         buffer_send%indx_offset(peer + 1) = buffer_send%indx_offset(peer) + send_count(peer)
         buffer_rec%indx_offset(peer + 1) = buffer_rec%indx_offset(peer) + recv_count(peer)
         IF (PRESENT(send_blocks)) THEN
            buffer_send%indx_offset(peer + 1) = buffer_send%indx_offset(peer) + send_blocks(peer)
            buffer_rec%indx_offset(peer + 1) = buffer_rec%indx_offset(peer) + recv_blocks(peer)
         END IF
      END DO
      ALLOCATE (buffer_send%msg(buffer_send%msg_offset(nprocs)), &
                buffer_rec%msg(buffer_rec%msg_offset(nprocs)), &
                buffer_send%indx(n_index, buffer_send%indx_offset(nprocs)), &
                buffer_rec%indx(n_index, buffer_rec%indx_offset(nprocs)))
   END SUBROUTINE allocate_matrix_buffers

! **************************************************************************************************
!> \brief Synchronously exchange flat matrix buffers using variable-count all-to-all.
!> \param para_env communicator covering all sending and receiving ranks
!> \param num_entries_rec numerical payload lengths received from each rank
!> \param num_entries_send numerical payload lengths sent to each rank
!> \param buffer_rec preallocated receive buffer
!> \param buffer_send populated send buffer
!> \param do_indx whether to exchange integer metadata (default: true)
!> \param do_msg whether to exchange numerical payloads (default: true)
! **************************************************************************************************
   SUBROUTINE exchange_matrix_buffers(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, &
                                      do_indx, do_msg)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: num_entries_rec, num_entries_send
      TYPE(matrix_buffer_type), INTENT(INOUT)            :: buffer_rec, buffer_send
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_indx, do_msg

      CHARACTER(LEN=*), PARAMETER :: routineN = 'exchange_matrix_buffers'

      INTEGER                                            :: handle, n_index, nprocs
      LOGICAL                                            :: my_do_indx, my_do_msg

      CALL timeset(routineN, handle)
      nprocs = para_env%num_pe
      my_do_indx = .TRUE.
      IF (PRESENT(do_indx)) my_do_indx = do_indx
      my_do_msg = .TRUE.
      IF (PRESENT(do_msg)) my_do_msg = do_msg

      IF (my_do_indx) THEN
         n_index = SIZE(buffer_send%indx, 1)
         CALL para_env%alltoall(buffer_send%indx, &
                                n_index*(buffer_send%indx_offset(1:) - buffer_send%indx_offset(:nprocs - 1)), &
                                n_index*buffer_send%indx_offset(:nprocs - 1), buffer_rec%indx, &
                                n_index*(buffer_rec%indx_offset(1:) - buffer_rec%indx_offset(:nprocs - 1)), &
                                n_index*buffer_rec%indx_offset(:nprocs - 1))
      END IF
      IF (my_do_msg) THEN
         CALL para_env%alltoall(buffer_send%msg, num_entries_send, buffer_send%msg_offset(:nprocs - 1), &
                                buffer_rec%msg, num_entries_rec, buffer_rec%msg_offset(:nprocs - 1))
      END IF
      CALL timestop(handle)
   END SUBROUTINE exchange_matrix_buffers

END MODULE matrix_communication
